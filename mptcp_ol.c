/*
 *	MPTCP Scheduler to reduce latency and jitter.
 *
 *	This scheduler sends all packets redundantly on all available subflows.
 *
 *	Initial Design & Implementation:
 *	Tobias Erbshaeusser <erbshauesser@dvs.tu-darmstadt.de>
 *	Alexander Froemmgen <froemmge@dvs.tu-darmstadt.de>
 *
 *	Initial corrections & modifications:
 *	Christian Pinedo <christian.pinedo@ehu.eus>
 *	Igor Lopez <igor.lopez@ehu.eus>
 *
 *	This program is free software; you can redistribute it and/or
 *      modify it under the terms of the GNU General Public License
 *      as published by the Free Software Foundation; either version
 *      2 of the License, or (at your option) any later version.
 */

#include <linux/module.h>
#include <net/mptcp.h>

#include <linux/random.h>

/* ytxing: yue */
typedef __u8 u8;
typedef __u16 u16;
typedef __u64 u64;
typedef int s32;

#define OLSCHED_INTERVALS_NUM 4
#define OLSCHED_INTERVALS_DURATION 1.5 /* OLSCHED_INTERVALS_DURATION * srtt_us */
#define OLSCHED_INTERVALS_TIMEOUT 10 /* OLSCHED_INTERVALS_TIMEOUT * srtt_us */
#define OLSCHED_INTERVAL_MIN_PACKETS 10

#define OLSCHED_SCALE 12
#define OLSCHED_UNIT (1 << OLSCHED_SCALE)
#define OLSCHED_MAX_RED_RATIO OLSCHED_UNIT
#define OLSCHED_MIN_RED_RATIO 0
#define OLSCHED_INIT_RED_RATIO 0
#define OLSCHED_PROBING_EPSILON (OLSCHED_UNIT * 0.01)

/* ytxing: for utility function calculation */
#define OLSCHED_ALPHA 0.99 /* (0, 1) */
#define OLSCHED_BETA 2
#define OLSCHED_GAMMA 2

/* Scale factor for rate in pkt/uSec unit to avoid truncation in bandwidth
 * estimation. The rate unit ~= (1500 bytes / 1 usec / 2^24) ~= 715 bps.
 * This handles bandwidths from 0.06pps (715bps) to 256Mpps (3Tbps) in a u32.
 * Since the minimum window is >=4 packets, the lower bound isn't
 * an issue. The upper bound isn't an issue with existing technologies.
 */
// #define BW_SCALE 24
// #define BW_UNIT (1 << BW_SCALE)



/* ytxing:	
 * ol_interval for each subflow.
 * This store some status at the beginning and the end of an interval,
 * in order to record the behaviour of this subflow and for utility 
 * function caculation.
 */
struct ol_interval {

	u8 index;
	u16	red_ratio;
	/* status for sending state */	
	u64	snd_time_begin;
	u64	snd_time_end;
	u32	pkts_out_begin;
	u32	pkts_out_end;
	u32	snd_seq_begin;
	u32	snd_seq_end;
	u32 delivered_begin; 	/*pkts*/
	u32 lost_begin; 		/*pkts*/

	/* status for receiving state */
	u64	rcv_begin_time;
	u64	rcv_end_time;
	u32 known_seq;			/* known send sequence (through (s)ack) */
	u32 delivered_end; 		/*pkts*/
	u32 lost_end; 			/*pkts*/

	s32 utility;

	u8	snd_timeout:1,
			snd_ended:1,
			rcv_timeout:1,
			rcv_ended:1,
			unusued:4;
};

enum ol_state {
	OL_PROBE, // sending data for probing 
	OL_WAIT,  // waiting for ack of sent data
	OL_MOVE,  // move toward one direction, stop when utility decreases
	OL_STAY,  // stay at current ratio
};

struct ol_global {
	u16	red_quota;	// ytxing: how many redundant pkts should this subflow send
	u16	new_quota;	// ytxing: how many new pkts should this subflow send
	u16	red_ratio;	// ytxing: (0, 128), multiply and >> OLSCHED_SCALE
	u16	suggested_red_ratio;	// ytxing: this ratio is suggested by current probing or decisions, and should be used to setup intervals lately

	/* intervals info */
	u8	snd_idx:4, 	// ytxing: the interval index currently sending pkts (< OLSCHED_INTERVALS_NUM)
		rcv_idx:4;	// ytxing: the interval index currently receiving pkts (< OLSCHED_INTERVALS_NUM)

	u8	waiting:1,	// ytxing: ture if all intervals finish sending but not stop receiving
		moving:1,	// ytxing: ture if olsched move toward one direction, stop when utility decreases
		unused:7;
	enum ol_state state;
};

/* Struct to store the data of a single subflow */
struct olsched_priv {
	/* The skb or NULL */
	struct sk_buff *skb;
	/* End sequence number of the skb. This number should be checked
	 * to be valid before the skb field is used
	 */
	u32 skb_end_seq;

	/* Structure to store status for each subflow in each interval */
	struct ol_interval *intervals_data;

	/* Some status for each subflow to follow currently */
	struct ol_global *global_data;
};

/* Struct to store the data of the control block */
struct olsched_cb {
	/* The next subflow where a skb should be sent or NULL */
	struct tcp_sock *next_subflow;
};

/* Returns the socket data from a given subflow socket */
static struct olsched_priv *olsched_get_priv(struct tcp_sock *tp)
{
	return (struct olsched_priv *)&tp->mptcp->mptcp_sched[0];
}

/* Returns the control block data from a given meta socket */
static struct olsched_cb *olsched_get_cb(struct tcp_sock *tp)
{
	return (struct olsched_cb *)&tp->mpcb->mptcp_sched[0];
}

static int olsched_get_active_valid_sks_num(struct sock *meta_sk)
{
	struct tcp_sock *meta_tp = tcp_sk(meta_sk);
	struct mptcp_cb *mpcb = meta_tp->mpcb;
	struct mptcp_tcp_sock *mptcp;
	int active_valid_sks = 0;
	int i = 0;

	mptcp_for_each_sub(mpcb, mptcp) {
		struct sock *sk = mptcp_to_sock(mptcp);
		i++;
		if (subflow_is_active((struct tcp_sock *)sk) &&
		    !mptcp_is_def_unavailable(sk))
			active_valid_sks++;
	}
		
	printk(KERN_DEBUG "ytxing: active_valid_sks: %d/%d\n", active_valid_sks, i);
	return active_valid_sks;
}


/* ytxing:	1. bandwidth = cwnd * mss / rtt 
 *			2. bandwidth = delivered_byte / time_interval
 *			in Bps
 */
static u64 olsched_get_bandwidth(struct sock *sk)
{
	struct tcp_sock *tp = tcp_sk(sk);
	u32 rtt_us, mss_now; 
	u64 bandwidth;

	if (tp->srtt_us) {		/* any RTT sample yet? */
		rtt_us = max(tp->srtt_us >> 3, 1U);
	} else {			 /* no RTT sample yet */
		rtt_us = USEC_PER_MSEC;	 /* use nominal default RTT ytxing maybe TODO */
	}
	
	mss_now = tcp_current_mss(sk);
	bandwidth = (u64)tp->snd_cwnd * mss_now * USEC_PER_SEC;
	do_div(bandwidth, rtt_us);
		
	printk(KERN_DEBUG "ytxing: tp: %p bandwidth%llu rtt:%u(us) cwnd:%u\n", tp, bandwidth, rtt_us, tp->snd_cwnd);
	return bandwidth; /* Bps */
}

/* ytxing:	
 * get the loss rate of the sk, make sure this interval is COMPLETED		
 * in (0, 1) << OLSCHED_SCALE
 */
static u64 olsched_get_loss_rate(struct sock *sk)
{
	struct tcp_sock *tp = tcp_sk(sk);
	struct olsched_priv *ol_p = olsched_get_priv(tp);
	struct ol_interval *interval = ol_p->intervals_data;
	u64 loss_rate;
	u32 delivered;

	/* nothing delivered*/
	if (interval->delivered_end == interval->delivered_begin)
		return (1 << OLSCHED_SCALE); 

	loss_rate = (interval->lost_end - interval->lost_begin) << OLSCHED_SCALE;
	delivered = interval->delivered_end - interval->delivered_begin;
	do_div(loss_rate, delivered);
	return loss_rate;
}


static bool olsched_use_subflow(struct sock *meta_sk,
				 int active_valid_sks,
				 struct tcp_sock *tp,
				 struct sk_buff *skb)
{
	if (!skb || !mptcp_is_available((struct sock *)tp, skb, false))
		return false;

	if (TCP_SKB_CB(skb)->path_mask != 0)
		return subflow_is_active(tp);

	if (TCP_SKB_CB(skb)->path_mask == 0) {
		if (active_valid_sks == -1)
			active_valid_sks = olsched_get_active_valid_sks_num(meta_sk);

		if (subflow_is_backup(tp) && active_valid_sks > 0)
			return false;
		else
			return true;
	}

	return false;
}

#define mptcp_entry_next_rcu(__mptcp)						\
	hlist_entry_safe(rcu_dereference_raw(hlist_next_rcu(			\
		&(__mptcp)->node)), struct mptcp_tcp_sock, node)

static void olsched_update_next_subflow(struct tcp_sock *tp,
					 struct olsched_cb *ol_cb)
{
	struct mptcp_tcp_sock *mptcp = mptcp_entry_next_rcu(tp->mptcp);

	if (mptcp)
		ol_cb->next_subflow = mptcp->tp;
	else
		ol_cb->next_subflow = NULL;
}

static struct sock *ol_get_available_subflow(struct sock *meta_sk,
					      struct sk_buff *skb,
					      bool zero_wnd_test)
{
	struct tcp_sock *meta_tp = tcp_sk(meta_sk);
	struct mptcp_cb *mpcb = meta_tp->mpcb;
	struct olsched_cb *ol_cb = olsched_get_cb(meta_tp);
	struct tcp_sock *first_tp = ol_cb->next_subflow, *tp;
	struct mptcp_tcp_sock *mptcp;
	int found = 0;

	/* Answer data_fin on same subflow */
	if (meta_sk->sk_shutdown & RCV_SHUTDOWN &&
	    skb && mptcp_is_data_fin(skb)) {
		mptcp_for_each_sub(mpcb, mptcp) {
			struct sock *sk = mptcp_to_sock(mptcp);

			if (tcp_sk(sk)->mptcp->path_index ==
				mpcb->dfin_path_index &&
			    mptcp_is_available(sk, skb, zero_wnd_test))
				return sk;
		}
	}

	if (!first_tp && !hlist_empty(&mpcb->conn_list)) {
		first_tp = hlist_entry_safe(rcu_dereference_raw(hlist_first_rcu(&mpcb->conn_list)),
					    struct mptcp_tcp_sock, node)->tp;
	}
	tp = first_tp;

	/* still NULL (no subflow in conn_list?) */
	if (!first_tp)
		return NULL;

	/* Search for a subflow to send it.
	 *
	 * We want to pick a subflow that is after 'first_tp' in the list of subflows.
	 * Thus, the first mptcp_for_each_sub()-loop tries to walk the list up
	 * to the subflow 'tp' and then checks whether any one of the remaining
	 * ones is eligible to send.
	 * The second mptcp_for_each-sub()-loop is then iterating from the
	 * beginning of the list up to 'first_tp'.
	 */
	mptcp_for_each_sub(mpcb, mptcp) {
		/* We go up to the subflow 'tp' and start from there */
		if (tp == mptcp->tp)
			found = 1;

		if (!found)
			continue;
		tp = mptcp->tp;

		if (mptcp_is_available((struct sock *)tp, skb,
				       zero_wnd_test)) {
			olsched_update_next_subflow(tp, ol_cb);
			return (struct sock *)tp;
		}
	}

	mptcp_for_each_sub(mpcb, mptcp) {
		tp = mptcp->tp;

		if (tp == first_tp)
			break;

		if (mptcp_is_available((struct sock *)tp, skb,
				       zero_wnd_test)) {
			olsched_update_next_subflow(tp, ol_cb);
			return (struct sock *)tp;
		}
	}

	/* No space */
	return NULL;
}

/* Corrects the stored skb pointers if they are invalid */
static void olsched_correct_skb_pointers(struct sock *meta_sk,
					  struct olsched_priv *ol_p)
{
	struct tcp_sock *meta_tp = tcp_sk(meta_sk);

	if (ol_p->skb &&
	    (!after(ol_p->skb_end_seq, meta_tp->snd_una) ||
	     after(ol_p->skb_end_seq, meta_tp->snd_nxt)))
		ol_p->skb = NULL;
}

/*	ytxing: get a bettersk that satisfies: 
 *			bestp->srtt < bettersk->srtt <= others->srtt,
 *			if bestp == NULL, then it simply returns the sk with min rtt.
 *			Should check if this sk can send a skb later. */
static struct sock *get_shorter_rtt_subflow_with_selector(struct mptcp_cb *mpcb, 
						struct sock *bestp, 
						bool (*selector)(const struct tcp_sock *))
{
	struct sock *bettersk = NULL;
	u32 better_srtt = 0xffffffff;
	u32 best_srtt = 0;
	struct mptcp_tcp_sock *mptcp;

	if (bestp){
		best_srtt = tcp_sk(bestp)->srtt_us;
	}

	mptcp_for_each_sub(mpcb, mptcp) {
		struct sock *sk = mptcp_to_sock(mptcp);
		struct tcp_sock *tp = tcp_sk(sk);

		/* First, we choose only the wanted sks */
		if (!(*selector)(tp))
			continue;

		if (mptcp_is_def_unavailable(sk))
			continue;

		if (tp->srtt_us <= best_srtt && bestp)
			continue;

		if (tp->srtt_us < better_srtt) {
			better_srtt = tp->srtt_us;
			bettersk = sk;
		}
	}

	return bettersk;
}

/* ytxing:	
 * Calculate the redundant quota of this subflow,
 * e.g., how many redundant should it send.
 * force means we go to another interval and refresh the quota accordingly.
 */
static void olsched_check_quota(struct tcp_sock *tp,
					  struct olsched_priv *ol_p, bool force)
{
	struct ol_global *global = ol_p->global_data;
	// still has quota, no need to calculate again
	if (!force && (ol_p->red_quota > 0 || ol_p->new_quota > 0))
		return;
	global->red_quota =  (tp->snd_cwnd * global->red_ratio) >> OLSCHED_SCALE;
	global->new_quota = tp->snd_cwnd - global->red_quota;
	printk(KERN_INFO "ytxing: tp %p red_quota %u new_quota %u ratio %u/255\n", tp, global->red_quota, global->new_quota, global->red_ratio);
}


/* Returns the next skb from the queue */
static struct sk_buff *olsched_next_skb_from_queue(struct sk_buff_head *queue,
						    struct sk_buff *previous,
						    struct sock *meta_sk)
{
	struct sk_buff *skb;

	if (!previous)
		return tcp_rtx_queue_head(meta_sk) ? : skb_peek(queue); 
		// skb_peek(queue) is the same as tcp_send_head(meta_sk)

	/* sk_data->skb stores the last scheduled packet for this subflow.
	 * If sk_data->skb was scheduled but not sent (e.g., due to nagle),
	 * we have to schedule it again.
	 *
	 * For the redundant scheduler, there are two cases:
	 * 1. sk_data->skb was not sent on another subflow:
	 *    we have to schedule it again to ensure that we do not
	 *    skip this packet.
	 * 2. sk_data->skb was already sent on another subflow:
	 *    with regard to the redundant semantic, we have to
	 *    schedule it again. However, we keep it simple and ignore it,
	 *    as it was already sent by another subflow.
	 *    This might be changed in the future.
	 *
	 * For case 1, send_head is equal previous, as only a single
	 * packet can be skipped.
	 */
	if (tcp_send_head(meta_sk) == previous)
		return tcp_send_head(meta_sk);// ytxing: this is the normal skb

	skb = skb_rb_next(previous);// ytxing: this is the redundant skb
	if (skb)
		return skb;

	return tcp_send_head(meta_sk);
}

/* ytxing:
 * calculate the utility of curr_sk
 * This function is only called once when one Scheduling Interval(SI) finishes.
 * When calling this function, there should be valid bandwidth 
 * and rtt sample for all subflows .
 */
#define FIXEDPT_BITS 64
#define FIXEDPT_WBITS 40
#define OMIT_STDINT
#include "fixedptc.h"
static u64 ol_calc_utility(struct sock *meta_sk, struct sock *curr_sk) {
	struct tcp_sock *meta_tp = tcp_sk(meta_sk), *curr_tp = tcp_sk(curr_sk);	
	struct mptcp_tcp_sock *mptcp;
	struct mptcp_cb *mpcb = meta_tp->mpcb;
	struct olsched_priv *curr_ol_p = olsched_get_priv(curr_tp);
	struct olsched_cb *ol_cb = olsched_get_cb(meta_tp);

	u64 others_bw, curr_bw, dRTT /* ((curr_rtt - minRTT) << OLSCHED_SCALE / minRTT) */;
	u32 ratio_t, minRTT = 0xFFFFFFFF;

	u64 bw_item, bw_item_f, ofo_item, loss_rate, loss_item, utility;

	/* ytxing:	in this loop, we 
	 *			1. calculate the ratio*bandwidth of other subflows 
	 *			2. find a minimal RTT
	 */
	mptcp_for_each_sub(mpcb, mptcp) {
		struct sock *sk = mptcp_to_sock(mptcp);		
		struct tcp_sock *tp = mptcp->tp;
		struct olsched_priv *ol_p;
		u64 bw_t;

		if (!subflow_is_active(tp) || mptcp_is_def_unavailable(sk)) /* we ignore unavailable subflows*/
			continue;

		if (tp->srtt_us < minRTT){
			minRTT = tp->srtt_us;
		}

		if (curr_sk == sk)
			continue;

		ol_p = olsched_get_priv(tp);
		ratio_t = ol_p->red_ratio;
		bw_t = olsched_get_bandwidth(sk); 
		others_bw = (1 - ratio_t) * bw_t; /* ytxing: (1 - ratio_t), be careful */
	}

	/* bandwidth item */
	curr_bw = olsched_get_bandwidth(curr_sk);
	bw_item = others_bw + curr_bw * curr_ol_p->red_ratio;
	printk(KERN_DEBUG "ytxing: tp:%p bw_item: %llu\n", curr_tp, bw_item);

	bw_item_f = fixedpt_fromint(bw_item);
	bw_item_f = fixedpt_pow(bw_item_f, fixedpt_rconst(OLSCHED_ALPHA));
	bw_item = (bw_item_f >> FIXEDPT_FBITS); /* OLSCHED_ALPHA power TODO (does it work?*/;
	printk(KERN_DEBUG "ytxing: tp:%p bw_item^%f: %llu\n", curr_tp, OLSCHED_ALPHA, bw_item);

	/* dRTT(ofo) item */
	if(curr_tp->srtt_us == minRTT) {
		dRTT = 0;
	}
	else {
		dRTT = (curr_tp->srtt_us - minRTT) << OLSCHED_SCALE;
		do_div(dRTT, minRTT);
		// dRTT = div_u64(dRTT, minRTT);
	}
	ofo_item = (bw_item * dRTT) >> OLSCHED_SCALE;

	/* loss item */
	loss_rate = olsched_get_loss_rate(curr_sk);
	loss_item = (bw_item * loss_rate) >> OLSCHED_SCALE; /* x loss_rate TODO*/

	/* finally, the utility function */
	utility = bw_item + OLSCHED_BETA * ofo_item + OLSCHED_GAMMA * loss_item;

	return utility;
}


static struct sk_buff *mptcp_ol_next_segment(struct sock *meta_sk,
					      int *reinject,
					      struct sock **subsk,
					      unsigned int *limit)
{
	struct tcp_sock *meta_tp = tcp_sk(meta_sk);
	struct mptcp_cb *mpcb = meta_tp->mpcb;
	struct olsched_cb *ol_cb = olsched_get_cb(meta_tp);
	struct tcp_sock *first_tp = ol_cb->next_subflow, *tp;
	struct mptcp_tcp_sock *mptcp;
	int active_valid_sks = -1;
	struct sk_buff *skb;
	int found = 0;

	/* As we set it, we have to reset it as well. */
	*limit = 0;

	if (skb_queue_empty(&mpcb->reinject_queue) &&
	    skb_queue_empty(&meta_sk->sk_write_queue) &&
	    tcp_rtx_queue_empty(meta_sk))
		/* Nothing to send */
		return NULL;

	/* First try reinjections */
	skb = skb_peek(&mpcb->reinject_queue);
	if (skb) {
		*subsk = get_available_subflow(meta_sk, skb, false);
		if (!*subsk)
			return NULL;
		*reinject = 1;
		return skb;
	}

	/* Then try indistinctly redundant and normal skbs */

	if (!first_tp && !hlist_empty(&mpcb->conn_list)) {
		first_tp = hlist_entry_safe(rcu_dereference_raw(hlist_first_rcu(&mpcb->conn_list)),
					    struct mptcp_tcp_sock, node)->tp;
	}

	/* still NULL (no subflow in conn_list?) */
	if (!first_tp)
		return NULL;

	tp = first_tp;

	*reinject = 0;
	active_valid_sks = olsched_get_active_valid_sks_num(meta_sk);

	/* We want to pick a subflow that is after 'first_tp' in the list of subflows.
	 * Thus, the first mptcp_for_each_sub()-loop tries to walk the list up
	 * to the subflow 'tp' and then checks whether any one of the remaining
	 * ones can send a segment.
	 * The second mptcp_for_each-sub()-loop is then iterating from the
	 * beginning of the list up to 'first_tp'.
	 */
	mptcp_for_each_sub(mpcb, mptcp) {
		struct olsched_priv *ol_p;

		if (tp == mptcp->tp)
			found = 1;

		if (!found)
			continue;

		tp = mptcp->tp;

		/* Correct the skb pointers of the current subflow */
		ol_p = olsched_get_priv(tp);
		olsched_correct_skb_pointers(meta_sk, ol_p);

		skb = olsched_next_skb_from_queue(&meta_sk->sk_write_queue,
						   ol_p->skb, meta_sk);
		if (skb && olsched_use_subflow(meta_sk, active_valid_sks, tp,
						skb)) {
			ol_p->skb = skb;
			ol_p->skb_end_seq = TCP_SKB_CB(skb)->end_seq;
			olsched_update_next_subflow(tp, ol_cb);
			*subsk = (struct sock *)tp;

			if (TCP_SKB_CB(skb)->path_mask)
				*reinject = -1;
			return skb;
		}
	}

	mptcp_for_each_sub(mpcb, mptcp) {
		struct olsched_priv *ol_p;

		tp = mptcp->tp;

		if (tp == first_tp)
			break;

		/* Correct the skb pointers of the current subflow */
		ol_p = olsched_get_priv(tp);
		olsched_correct_skb_pointers(meta_sk, ol_p);

		skb = olsched_next_skb_from_queue(&meta_sk->sk_write_queue,
						   ol_p->skb, meta_sk);
		if (skb && olsched_use_subflow(meta_sk, active_valid_sks, tp,
						skb)) {
			ol_p->skb = skb;
			ol_p->skb_end_seq = TCP_SKB_CB(skb)->end_seq;
			olsched_update_next_subflow(tp, ol_cb);
			*subsk = (struct sock *)tp;

			if (TCP_SKB_CB(skb)->path_mask)
				*reinject = -1;
			return skb;
		}
	}

	/* Nothing to send */
	return NULL;
}

// ytxing: 	1. find a subflow with minimal RTT;
//			2. get a redundant or normal skb for the subflow
static struct sk_buff *mptcp_ol_next_segment_rtt(struct sock *meta_sk,
					      int *reinject,
					      struct sock **subsk,
					      unsigned int *limit)
{
	struct tcp_sock *meta_tp = tcp_sk(meta_sk);
	struct mptcp_cb *mpcb = meta_tp->mpcb;
	struct olsched_cb *ol_cb = olsched_get_cb(meta_tp);
	struct olsched_priv *ol_p;
	struct tcp_sock *tp;
	// struct mptcp_tcp_sock *mptcp;
	int active_valid_sks = -1;
	struct sk_buff *skb;
	struct sock *sk, *unavailable_sk;
	// int found = 0;
	int i = 0;

	/* As we set it, we have to reset it as well. */
	*limit = 0;

	if (skb_queue_empty(&mpcb->reinject_queue) &&
	    skb_queue_empty(&meta_sk->sk_write_queue) &&
	    tcp_rtx_queue_empty(meta_sk))
		/* Nothing to send */
		return NULL;

	/* First try reinjections */
	skb = skb_peek(&mpcb->reinject_queue);
	if (skb) {
		*subsk = get_available_subflow(meta_sk, skb, false);
		if (!*subsk)
			return NULL;
		*reinject = 1;
		return skb;
	}

	/* Then try indistinctly redundant and normal skbs */

	*reinject = 0;
	active_valid_sks = olsched_get_active_valid_sks_num(meta_sk);

	unavailable_sk = NULL;
	/* ytxing: in this loop we find a avaliable sk with rtt as short as possible*/
	for(i = 0; i < active_valid_sks; i++){
		printk(KERN_DEBUG "ytxing: i: %d/%d\n", i+1, active_valid_sks);
		sk = get_shorter_rtt_subflow_with_selector(mpcb, unavailable_sk ,&subflow_is_active);
		if (!sk)
			break;
		tp = tcp_sk(sk);
		printk(KERN_DEBUG "ytxing: better tp: %p srtt: %u\n", tp, (tp->srtt_us >> 3));
		ol_p = olsched_get_priv(tp);
		olsched_correct_skb_pointers(meta_sk, ol_p);

		/* process ol model here? */
		olsched_check_quota(tp, ol_p); // ytxing
		if (ol_p->red_quota >= 1){ // Executing redundant routines
			skb = olsched_next_skb_from_queue(&meta_sk->sk_write_queue,
								ol_p->skb, meta_sk);
		
			if (skb && olsched_use_subflow(meta_sk, active_valid_sks, tp,
							skb)) {
				ol_p->skb = skb;
				ol_p->skb_end_seq = TCP_SKB_CB(skb)->end_seq;
				olsched_update_next_subflow(tp, ol_cb);
				*subsk = (struct sock *)tp;
				// printk(KERN_DEBUG "ytxing skb_len %u (mss?)\n", skb->len);
				ol_p->red_quota -= 1;
				if (TCP_SKB_CB(skb)->path_mask){
					*reinject = -1;
					// printk(KERN_DEBUG "ytxing: a red pkt\n");
					}
				// printk(KERN_DEBUG "ytxing: send pkt 1\n");
				return skb;
			}
		}
		else { // send a new skb, no need redundant
			skb = tcp_send_head(meta_sk);
			if (skb && olsched_use_subflow(meta_sk, active_valid_sks, tp,
							skb)) {
				ol_p->skb = skb;
				ol_p->skb_end_seq = TCP_SKB_CB(skb)->end_seq;
				olsched_update_next_subflow(tp, ol_cb);
				*subsk = (struct sock *)tp;
				// printk(KERN_DEBUG "ytxing skb_len %u (mss?)\n", skb->len);
				ol_p->new_quota -= 1;
				if (TCP_SKB_CB(skb)->path_mask)
					*reinject = -1;
				// printk(KERN_DEBUG "ytxing: send pkt 2\n");	
				return skb;
			}
		}

		// ytxing: 	now this sk is unable to send (maybe cwnd limited)
		//			to a sk with longer rtt.
		unavailable_sk = sk;
	}
	// printk(KERN_DEBUG "ytxing: tp: %p Nothing to send\n",  tp);

	/* Nothing to send */
	return NULL;
}

/******************
 * Intervals init *
 * ****************/

 /* Set the target rates of all intervals and reset statistics. */
static void ol_setup_intervals_probing(struct tcp_sock *tp)
{	
	struct olsched_priv *ol_p;
	ol_p = olsched_get_priv(tp);
	struct ol_global *global = ol_p->global_data;
	struct ol_interval *interval;
	

	u16 red_ratio_hi, red_ratio_lo;
	char rand;
	int i;

	get_random_bytes(&rand, 1);

	/* change the global ratio for probing */
	red_ratio_hi = min(global->suggested_red_ratio + OLSCHED_PROBING_EPSILON, OLSCHED_MAX_RED_RATIO);
	red_ratio_lo = max(global->suggested_red_ratio - OLSCHED_PROBING_EPSILON, OLSCHED_MIN_RED_RATIO);

	for (i = 0; i < OLSCHED_INTERVALS_NUM; i += 2) {
		if ((rand >> (i / 2)) & 1) {
			ol_p->intervals_data[i].red_ratio = red_ratio_lo;
			ol_p->intervals_data[i + 1].red_ratio = red_ratio_hi;
		} else {
			ol_p->intervals_data[i].red_ratio = red_ratio_hi;
			ol_p->intervals_data[i + 1].red_ratio = red_ratio_lo;
		}

		ol_p->intervals_data[i].index = i;
		ol_p->intervals_data[i + 1].index = i + 1;

	}

	global->snd_idx = 0;
	global->rcv_idx = 0;
	global->waiting = false;
	global->state = OL_PROBE;
}

/* Reset statistics and set the target rate for just one monitor interval */
static void ol_setup_intervals_moving(struct tcp_sock *tp)
{
	struct olsched_priv *ol_p;
	ol_p = olsched_get_priv(tp);
	struct ol_global *global = ol_p->global_data;

	ol_p->intervals_data[0].red_ratio = ol_p->global_data->suggested_red_ratio;
	global->snd_idx = 0;
	global->rcv_idx = 0;
	global->waiting = false;
	global->state = OL_MOVE;
}

/* Set the red_ratio based on the currently-sending interval
 * and update some interval states.
 */
static void start_current_send_interval(struct tcp_sock *tp)
{
	struct olsched_priv *ol_p;
	ol_p = olsched_get_priv(tp);
	struct ol_global *global = ol_p->global_data;
	struct ol_interval *interval = ol_p->intervals_data[global->snd_idx];

	interval->snd_time_begin = tp->tcp_mstamp;
	interval->snd_seq_begin = tp->snd_next;
	interval->pkts_out_begin = tp->data_segs_out; /* maybe? */
	interval->known_seq = tp->snd_una; /* init known_seq as the next unacked seq, should be less than snd_next*/

	/* when a interval is active, the global ratio is set to the interval's ratio */
	global->red_ratio = interval->red_ratio;

}


/* Have we sent all the data we need to for this interval? Must have at least
 * the minimum number of packets and should have sent 1 RTT worth of data.
 */
bool ol_current_send_interval_ended(struct ol_interval *interval, struct tcp_sock *tp)
{
	u32 packets_sent = tp->data_segs_out - interval->pkts_out_begin;
	u32 srtt_us = tp->srtt_us >> 3;
	bool timeout;

	/* not enough sending duration */
	if (tp->tcp_mstamp - interval->snd_time_begin < srtt_us * OLSCHED_INTERVALS_DURATION) 
		return false;
		
	timeout = (tp->tcp_mstamp - interval->snd_time_begin) > srtt_us * OLSCHED_INTERVALS_TIMEOUT;
	/* not enough packets out and not timeout */
	if (packets_sent < OLSCHED_INTERVAL_MIN_PACKETS && !timeout)
		return false;

	/* end the sending state this interval */
	interval->pkts_out_end = tp->data_segs_out; /* sure? maybe for non-sack */
	interval->snd_seq_end = tp->snd_nxt;
	interval->snd_time_end = tp->tcp_mstamp;
	interval->snd_timeout = timeout;
	interval->snd_ended = true;
	return true;
}

/* Have we accounted for (acked or lost) enough of the packets that we sent to
 * calculate summary statistics?
 */
bool recive_interval_ended(struct ol_interval *interval, struct tcp_sock *tp)
{
	if(interval->rcv_ended)
		return true;

	if(interval->snd_ended && !before(interval->known_seq, interval->snd_seq_end)){
		interval->rcv_ended = true;
		return true;
	}

	bool timeout = (tp->tcp_mstamp - interval->rcv_begin_time) > srtt_us * OLSCHED_INTERVALS_TIMEOUT;
	if(timeout){
		interval->rcv_timeout = timeout;
		interval->rcv_ended = true;
		if(interval->snd_ended == false)
			printk(KERN_DEBUG "ytxing: BUG tp: %p interval:%d rcv_timeout before snd_ended\n", tp, interval->index);
		return true;
	}

	return false;
}

bool all_recive_interval_ended(struct ol_global *global, struct tcp_sock *tp)
{
	int i;
	struct olsched_priv *ol_p;
	ol_p = olsched_get_priv(tp);
	for(i = 0; i < OLSCHED_INTERVALS_NUM < i++){
		if(!recive_interval_ended(ol_p->intervals_data[i], tp)){
			return false;
		}
	}
	return true;
}


/* ytxing:
 * change the state if needed, and process to another interval in OL_PROBE.
 */
static void interval_process(struct ol_global *global, struct tcp_sock *tp)
{
	global->snd_idx++;
	// if (global->snd_idx == OLSCHED_INTERVALS_NUM || global->moving) {
	// 	global->waiting = true;
	// }
	if ((global->snd_idx == OLSCHED_INTERVALS_NUM && global->state == OL_PROBE) 	/* all intervals stop sending in OL_PROBE */
		|| (global->snd_idx == 1 && global->state == OL_MOVE) 						/* the only one interval stop sending in OL_MOVE */
		|| (global->snd_idx == 1 && global->state == OL_STAY) {						/* the only one interval stop sending in OL_STAY */

		global->state = OL_WAIT;
		return;
	}

	start_current_send_interval(tp);
}

static void update_interval_info_snd(struct ol_interval *interval, struct tcp_sock *tp)
{
	interval->snd_seq_end = tp->snd_nxt;
}

/* ytxing:
 * update_interval_info_rcv collect infomation from sack blocks
 * in order to track the pkts sent in a certain interval.
 * When the interval sees a seq (known_seq) larger than the snd_seq_end, and the sending has stopped,
 * mark the time stamp.
 */
static void update_interval_info_rcv(struct ol_global *global, struct tcp_sock *tp)
{
	struct tcp_sock *tp = tcp_sk(sk);
	struct olsched_priv *ol_p;
	ol_p = olsched_get_priv(tp);
	int i,j;
	struct tcp_sack_block sack_cache[4];

	//sort received sacks according to sequence in increasing order
	if (tcp_is_sack(tp) && tp->sacked_out) {
		memcpy(sack_cache, tp->recv_sack_cache, sizeof(sack_cache));
		for (i = 0; i < 4; i++) {
			for (j = i+1; j < 4; j++) {
				if (after(sack_cache[i].start_seq, sack_cache[j].start_seq)) {
					u32 tmp = sack_cache[i].start_seq;
					sack_cache[i].start_seq = sack_cache[j].start_seq;
					sack_cache[j].start_seq = tmp;

					tmp = sack_cache[i].end_seq;
					sack_cache[i].end_seq = sack_cache[j].end_seq;
					sack_cache[j].end_seq = tmp;
				}
			}
		}
	}
	
	//for all active intervals check if cumulative acks changed the last known seq, or if the sacks did
	for (i = 0; i < OLSCHED_INTERVALS_NUM; i++) {
		struct ol_interval *interval = ol_p->intervals_data[i];
		u32 current_interval_known_seq = interval->known_seq;

		//set the last known sequence to the last cumulative ack if it is better than the last known seq
		if (after(tp->snd_una, interval->known_seq)) {
			interval->known_seq = tp->snd_una;
		}

		//there are sacks
		if (tcp_is_sack(tp) && tp->sacked_out) {
			for (j = 0; j < 4; j++) {
				//if the sack doesn't bring any new information, check the next one
				/* ytxing: then interval->rcv_ended should be true */
				if (!before(interval->known_seq, interval->end_seq)) {
					continue;
				}

				//mark the hole as lost bytes in this monitor interval
				if (sack_cache[j].start_seq != 0 && sack_cache[j].end_seq != 0) {
					if (before(interval->known_seq, sack_cache[j].start_seq)) {
						s32 lost = min(sack_cache[j].start_seq - interval->known_seq, interval->snd_seq_end - interval->known_seq);
					}
					//update the last known seq if it was changed
					if (after(sack_cache[j].end_seq, interval->known_seq)) {
						interval->known_seq = sack_cache[j].end_seq;
					}
				}
			}
		}
		/* ytxing: if this interval just start receiving pkts for the first time */
        if (before(current_interval_known_seq, interval->start_seq) && !before(interval->known_seq, interval->start_seq)) {
            interval->rcv_begin_time = tcp_sk(sk)->tcp_mstamp;
        }
		/* ytxing: if the interval has seen all sent packets and finished sending */
        if (interval->snd_ended && before(current_interval_known_seq, interval->snd_seq_end) && !before(interval->known_seq, interval->snd_seq_end)) {
            interval->rcv_end_time = tcp_sk(sk)->tcp_mstamp;
			interval->rcv_ended = true;
			// interval->delivered_end = tp->delivered; /* sure? delivered is not reliable with sack since it might contain some pkts of other intervals*/
        }

	}

}
/* Updates the OL model */
static void ol_process(struct sock *meta_sk, struct sock *sk)
{

	struct tcp_sock *meta_tp = tcp_sk(meta_sk), *tp = tcp_sk(sk);
	struct olsched_cb *ol_cb;
	ol_cb = olsched_get_cb(meta_tp);

	struct olsched_priv *ol_p;
	ol_p = olsched_get_priv(tp);
	struct ol_global *global = ol_p->global_data;
	struct ol_interval *interval;

	/* maybe some checking? */
	// if (!ol_valid(ol_p))
	// 	return;

	/* handle interval in sending state */
	if (!global->state == OL_WAIT) { 
		interval = &ol_p->intervals_data[global->snd_idx];
		update_interval_info_snd(interval, tp);
		if (ol_current_send_interval_ended(interval, tp)) {
			interval_process(global, tp);
		}
	}

	/* handle interval in receiving state */ //TODO 这里需要实时记录一些必要的数据在interval中
	update_interval_info_rcv(global, tp);
	/* check if all intervals finish receiving */
	if (all_recive_interval_ended(global, tp)){
		/* decide new suggested ratio and ol_state TODO */
		/* setup intervals again TODO */
	}
}

static void olsched_release(struct sock *sk)
{
	struct tcp_sock *tp = tcp_sk(sk);
	struct olsched_cb *ol_cb = olsched_get_cb(tp);

	/* Check if the next subflow would be the released one. If yes correct
	 * the pointer
	 */
	if (ol_cb->next_subflow == tp)
		olsched_update_next_subflow(tp, ol_cb);
}

static void olsched_init(struct sock *sk)
{
	struct olsched_priv *ol_p = olsched_get_priv(tcp_sk(sk));
	struct olsched_cb *ol_cb = olsched_get_cb(tcp_sk(mptcp_meta_sk(sk)));

	ol_p->intervals_data = kzalloc(sizeof(struct ol_interval) * OL_INTERVALS * 2,
				 GFP_KERNEL);
	ol_p->global_data = kzalloc(sizeof(struct ol_global),
				 GFP_KERNEL);
	ol_p->global_data->red_ratio = OLSCHED_UNIT * OLSCHED_INIT_RED_RATIO;
	ol_p->global_data->red_quota = 0;
	ol_p->global_data->new_quota = 0;

	/* first try to move from 0 -> 1, for throughput TODO */

	printk(KERN_DEBUG "ytxing: olsched_init tp %p\n", tcp_sk(sk));
	if(!tcp_is_sack)
		printk(KERN_DEBUG "ytxing: olsched_init tp %p is not sack, maybe bad.\n", tcp_sk(sk));
}

static struct mptcp_sched_ops mptcp_sched_ol = {
	.get_subflow = ol_get_available_subflow,
	.next_segment = mptcp_ol_next_segment_rtt,
	.init = olsched_init,
	.release = olsched_release,
	.name = "ol",
	.owner = THIS_MODULE,
};

static int __init ol_register(void)
{
	printk(KERN_DEBUG "ytxing: olsched_priv size :%lu (< %u)\n", sizeof(struct olsched_priv), MPTCP_SCHED_SIZE);
	printk(KERN_DEBUG "ytxing: olsched_cb size :%lu (< %u)\n", sizeof(struct olsched_cb), MPTCP_SCHED_DATA_SIZE);
	BUILD_BUG_ON(sizeof(struct olsched_priv) > MPTCP_SCHED_SIZE);
	BUILD_BUG_ON(sizeof(struct olsched_cb) > MPTCP_SCHED_DATA_SIZE);
	if (mptcp_register_scheduler(&mptcp_sched_ol))
		return -1;
	return 0;
}

static void ol_unregister(void)
{
	mptcp_unregister_scheduler(&mptcp_sched_ol);
}

module_init(ol_register);
module_exit(ol_unregister);

MODULE_AUTHOR("Yitao Xing");
MODULE_LICENSE("GPL");
MODULE_DESCRIPTION("OL MPTCP");
MODULE_VERSION("0.01");
