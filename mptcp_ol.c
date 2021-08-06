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
typedef __s8 s8;
typedef __s32 s32;
typedef __s64 s64;

#define OLSCHED_INTERVALS_NUM 1
#define OLSCHED_INTERVALS_MIN_DURATION 100 * USEC_PER_MSEC /* minimal duration(us) */
#define OLSCHED_INTERVALS_TIMEOUT OLSCHED_INTERVALS_MIN_DURATION * 3
#define OLSCHED_INTERVAL_MIN_PACKETS 30

#define OLSCHED_ARMS_NUM 4

#define OLSCHED_SCALE 13
#define OLSCHED_UNIT (1 << OLSCHED_SCALE)
#define OLSCHED_MAX_RED_RATIO OLSCHED_UNIT
#define OLSCHED_MIN_RED_RATIO 1
#define OLSCHED_INIT_RED_RATIO (1 * OLSCHED_UNIT)
#define OLSCHED_PROBING_EPSILON (OLSCHED_UNIT >> (OLSCHED_SCALE - 10))
#define OLSCHED_GRAD_NOISE (OLSCHED_UNIT >> (OLSCHED_SCALE - 3))
#define OLSCHED_GRAD_STEP 1/1000
#define OLSCHED_GRAD_MIN_STEP (OLSCHED_UNIT >> 4)
#define OLSCHED_START_STEP (OLSCHED_UNIT >> 3)


/* ytxing: for utility function calculation */
#define OLSCHED_GAMMA_MAB 99 /* div by OLSCHED_GAMMA_MAB_BASE */ 
#define OLSCHED_GAMMA_MAB_BASE 100

static const u16 ol_arm_to_red_ratio[OLSCHED_ARMS_NUM] = {
	OLSCHED_UNIT * 1,
	OLSCHED_UNIT * 2 / 3,
	OLSCHED_UNIT * 1 / 3,
	OLSCHED_MIN_RED_RATIO
};

/* Scale factor for rate in pkt/uSec unit to avoid truncation in bandwidth
 * estimation. The rate unit ~= (1500 bytes / 1 usec / 2^24) ~= 715 bps.
 * This handles bandwidths from 0.06pps (715bps) to 256Mpps (3Tbps) in a u32.
 * Since the minimum window is >=4 packets, the lower bound isn't
 * an issue. The upper bound isn't an issue with existing technologies.
 */
// #define BW_SCALE 24
// #define BW_UNIT (1 << BW_SCALE)

enum OL_MOVING_STATE {
	OL_RATIO_EQUAL, 
	OL_RATIO_UP,  
	OL_RATIO_DOWN,  
};

/* ytxing:	
 * ol_interval for each subflow.
 * This store some status at the beginning and the end of an interval,
 * in order to record the behaviour of this subflow and for utility 
 * function caculation.
 */
struct ol_interval_MAB {

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
	u64	rcv_time_begin;
	u64	rcv_time_end;
	u32 known_seq;			/* known send sequence (through (s)ack) */
	u32 delivered_end; 		/*pkts*/
	u32 lost_end; 			/*pkts*/
	u32 lost_bytes;

	u8 previous_arm_idx;
	u8 current_arm_idx;
	u16 arm_weight[OLSCHED_ARMS_NUM];
	u16 arm_probability[OLSCHED_ARMS_NUM];

	u32	interval_duration;

	u8	snd_ended:1,
		rcv_ended:1,
		invalid_utility:1,
		init:1,
		unusued:4;
};

enum OL_GLOBAL_STATE {
	OL_PROBE, // sending data for probing 
	OL_MOVE,  // move toward one direction, stop when utility decreases
	OL_START,  // stay at current ratio
};


struct ol_global_MAB {
	u16	red_quota;	// ytxing: how many redundant pkts should this subflow send
	u16	new_quota;	// ytxing: how many new pkts should this subflow send
	u16	red_ratio;	// ytxing: the ratio to calculate current quota, the same as the red_ratio current interval >> OLSCHED_SCALE
	

	struct ol_interval_MAB *intervals_data; // ytxing: for meta_tp infomation

	/* intervals info */
	u8	snd_idx:4, 	// ytxing: the interval index currently sending pkts (< OLSCHED_INTERVALS_NUM)
		rcv_idx:4;	// ytxing: the interval index currently receiving pkts (< OLSCHED_INTERVALS_NUM)

	u8	waiting:1,	// ytxing: ture if all intervals finish sending but not stop receiving
		moving:1,	// ytxing: ture if olsched move toward one direction, stop when utility decreases
		first_move:1,
		init:1,
		unused:4; 

	enum OL_GLOBAL_STATE state;	

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
	struct ol_interval_MAB *intervals_data;

	/* Some status for each subflow to follow currently */
	struct ol_global_MAB *global_data;
};

/* Struct to store the data of the control block */
struct olsched_cb {
	/* The next subflow where a skb should be sent or NULL */
	struct tcp_sock *next_subflow;
	struct ol_interval_MAB *meta_interval;
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

/* get x = number * OLSCHED_UNIT, return (e^number)*OLSCHED_UNIT */
static u32 ol_exp(u32 x)
{
	s64 temp = OLSCHED_UNIT;
	s64 e = OLSCHED_UNIT;
	int i;

	for (i = 1; temp != 0; i++) {
		temp *= x;
		temp /= i;
		temp /= OLSCHED_UNIT;
		e += temp;
	}
	return e;
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
		
	// printk(KERN_DEBUG "ytxing: active_valid_sks: %d/%d\n", active_valid_sks, i);
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
		return 0;
	}
	
	mss_now = tcp_current_mss(sk);
	bandwidth = (u64)tp->snd_cwnd * mss_now * USEC_PER_SEC;
	do_div(bandwidth, rtt_us);
		
	// printk(KERN_DEBUG "ytxing: tp:%p bandwidth%llu rtt:%u(us) cwnd:%u\n", tp, bandwidth, rtt_us, tp->snd_cwnd);
	return bandwidth; /* Bps */
}

/* ytxing:	1. bandwidth = cwnd * mss / rtt 
 *			2. bandwidth = delivered_byte / time_interval
 *			in Bps
 */
static u64 olsched_get_bandwidth_interval(struct ol_interval_MAB *interval, struct tcp_sock *tp)
{
	u64 bandwidth;
	u32 mss = tp->mss_cache;
	u64 byte_rcv = (interval->delivered_end - interval->delivered_begin) * mss;
	u64 duration = interval->rcv_time_end - interval->rcv_time_begin; 
	if (duration <= 0 || byte_rcv <= 0)
		return 0;

	bandwidth = byte_rcv * USEC_PER_SEC;
	do_div(bandwidth, duration); /* delivery rate, actually */
	return bandwidth; /* Bps */
}

/* ytxing:	
 * get the loss rate of the sk, make sure this interval is COMPLETED		
 * in (0, 1) << OLSCHED_SCALE
 */
static u64 olsched_get_loss_rate(struct ol_interval_MAB *interval)
{
	u64 loss_rate;
	u32 delivered;

	/* nothing delivered */
	if (interval->delivered_end == interval->delivered_begin)
		return (1 << OLSCHED_SCALE); 

	loss_rate = (interval->lost_end - interval->lost_begin) << OLSCHED_SCALE;
	delivered = interval->delivered_end - interval->delivered_begin;

	do_div(loss_rate, delivered);
	return loss_rate;
}

/* tweak the suggested ratio and setup four intervals */
static void ol_setup_intervals_MAB(struct tcp_sock *tp)
{	
	struct olsched_priv *ol_p = olsched_get_priv(tp);
	struct olsched_cb *ol_cb = olsched_get_cb(tp);
	struct ol_global_MAB *global = ol_p->global_data;
	struct ol_interval_MAB *interval, *meta_interval = ol_cb->meta_interval[0];
	
	

	global->snd_idx = 0;
	global->rcv_idx = 0;
	global->waiting = false;
	global->state = OL_PROBE;
}


/* tweak the suggested ratio and setup four intervals */
static void ol_setup_intervals_probing(struct tcp_sock *tp, int arm_idx)
{	
	struct olsched_priv *ol_p;
	struct ol_global_MAB *global;


	for (i = 0; i < OLSCHED_INTERVALS_NUM; i += 2) {
		if ((rand >> (i / 2)) & 1) {
			ol_p->intervals_data[i].red_ratio = red_ratio_lo;
			ol_p->intervals_data[i + 1].red_ratio = red_ratio_hi;
		} else {
			ol_p->intervals_data[i].red_ratio = red_ratio_hi;
			ol_p->intervals_data[i + 1].red_ratio = red_ratio_lo;
		}

	}

	global->snd_idx = 0;
	global->rcv_idx = 0;
	global->waiting = false;
	global->state = OL_PROBE;
}

/* follow the suggested ratio and setup one interval */
static void ol_setup_intervals_moving(struct tcp_sock *tp)
{
	struct olsched_priv *ol_p;
	struct ol_global *global;
	ol_p = olsched_get_priv(tp);
	global = ol_p->global_data;

	ol_p->intervals_data[0].red_ratio = ol_p->global_data->suggested_red_ratio;
	if(global->state == OL_PROBE) /* when switching from OL_PROBE to OL_MOVE, choose a higher utility for comparison */
		ol_p->intervals_data[0].previous_utility = max(ol_p->intervals_data[0].utility, ol_p->intervals_data[1].utility); 
	else /* when stay in OL_MOVE */
		ol_p->intervals_data[0].previous_utility = ol_p->intervals_data[0].utility;


	global->snd_idx = 0;
	global->rcv_idx = 0;
	global->waiting = false;
	global->state = OL_MOVE;
}

/* follow the suggested ratio and setup one interval */
static void ol_setup_intervals_starting(struct tcp_sock *tp)
{
	struct olsched_priv *ol_p;
	struct ol_global_MAB *global;
	ol_p = olsched_get_priv(tp);
	global = ol_p->global_data;

	ol_p->intervals_data[0].red_ratio = ol_p->global_data->suggested_red_ratio;
	ol_p->intervals_data[0].previous_utility = ol_p->intervals_data[0].utility; /* for comparison */

	global->snd_idx = 0;
	global->rcv_idx = 0;
	global->waiting = false;
	global->state = OL_START;
}

static u16 ol_suggest_ratio_probing(struct tcp_sock *tp)
{
	struct olsched_priv *ol_p;
	struct ol_global *global;
	bool run_1_res, run_2_res, did_agree;

	ol_p = olsched_get_priv(tp);
	global = ol_p->global_data;

	run_1_res = ol_p->intervals_data[0].utility > ol_p->intervals_data[1].utility;
	run_2_res = ol_p->intervals_data[2].utility > ol_p->intervals_data[3].utility;

	/* did_agree: was the 2 sets of intervals with the same result */
	did_agree = !((run_1_res == run_2_res) ^
				  (ol_p->intervals_data[0].red_ratio == ol_p->intervals_data[2].red_ratio)); 

	if (did_agree) {
		/* return the suggested ratio later */
		printk(KERN_INFO "ytxing: tp:%p OL_PROBE agreement ACHIEVED\n", tp);
		return run_2_res ? ol_p->intervals_data[2].red_ratio : ol_p->intervals_data[3].red_ratio;
	} else {
		/* no agreed result, return the same suggested ratio */
		printk(KERN_INFO "ytxing: tp:%p OL_PROBE agreement NOT achieved\n", tp);
		return global->suggested_red_ratio;
	}
}

/* Choose up/down rate changes based on utility gradient */
static u32 ol_suggest_ratio_moving(struct tcp_sock *tp)
{
	struct sock *sk = (struct sock *)tp;
	struct olsched_priv *ol_p = olsched_get_priv(tp) ;
	struct ol_global *global;
	struct ol_interval *interval = &ol_p->intervals_data[0];
	s64 utility, previous_utility, grad, grad_sign;
	u64 step_abs, grad_abs;
	s64 new_ratio, step;
	
	ol_p = olsched_get_priv(tp);
	global = ol_p->global_data;

	previous_utility = interval->utility;
	global->util_func(interval, sk);
	utility = interval->utility;


	grad = ol_calc_util_grad(global->suggested_red_ratio, utility, global->previous_suggested_red_ratio, previous_utility);
	printk(KERN_INFO "ytxing: tp:%p utility:%lld previous_utility:%lld\n", tp, utility, previous_utility);
	printk(KERN_INFO "ytxing: tp:%p global->suggested_red_ratio:%u global->previous_suggested_red_ratio:%u\n", tp, global->suggested_red_ratio >> 3, global->previous_suggested_red_ratio >> 3);
	if (grad == 0){	
		printk(KERN_INFO "ytxing: tp:%p OL_MOVE grad %lld stay\n", tp, grad);
		return global->suggested_red_ratio;
	}

	grad_abs = abs(grad);
	grad_sign = grad >= 0 ? 1 : -1; 
	step_abs = (grad_abs * OLSCHED_GRAD_STEP);  /* gradient ascent */
	step_abs = step_abs >> OLSCHED_SCALE;

	/* We need our step size to be large enough that we can compute the gradient
	 * with low noise
	 */
	step_abs = max_t(u64, step_abs, OLSCHED_GRAD_MIN_STEP);
	step = step_abs * grad_sign;

	new_ratio = max_t(s64, global->suggested_red_ratio + step, OLSCHED_MIN_RED_RATIO);
	new_ratio = min_t(s64, new_ratio, OLSCHED_MAX_RED_RATIO);

	printk(KERN_INFO "ytxing: tp:%p OL_MOVE grad %lld step %lld\n", tp, grad / 8, step / 8);

	return new_ratio;
}

/* Have we sent all the data we need to for this interval? Must have at least
 * the minimum number of packets and should have sent 1 RTT worth of data.
 */
bool ol_current_send_interval_ended(struct ol_interval_MAB *interval, struct tcp_sock *tp)
{
	u32 packets_sent = tp->data_segs_out - interval->pkts_out_begin;
	u32 interval_duration = interval->interval_duration;


	/* not enough sending duration */
	if (tp->tcp_mstamp - interval->snd_time_begin < interval_duration) 
		return false;
		
	// /* not enough packets out and not timeout */
	// if (packets_sent < OLSCHED_INTERVAL_MIN_PACKETS && !timeout)
	// 	return false;

	/* end the sending state this interval */
	interval->pkts_out_end = tp->data_segs_out; /* sure? maybe for non-sack */
	interval->snd_seq_end = tp->snd_nxt;
	interval->snd_time_end = tp->tcp_mstamp;
	interval->snd_ended = true;
	printk(KERN_INFO "ytxing: tp:%p idx:%d SND END packets_sent:%u\n", tp, interval->index, packets_sent);
	printk(KERN_INFO "ytxing: tp:%p idx:%d SND END actual_duration:%llu interval_duration:%u\n", tp, interval->index, interval->snd_time_end - interval->snd_time_begin, interval_duration);
	return true;
}

/* all send intervals (of all subflows) ended? */
bool ol_all_subflow_send_interval_ended(struct tcp_sock *meta_tp)
{
	struct mptcp_tcp_sock *mptcp;
	// struct ol_global *global = olsched_get_cb(meta_tp);
	struct mptcp_cb *mpcb = meta_tp->mpcb;
	bool all_ended = true;

	mptcp_for_each_sub(mpcb, mptcp) {
		struct sock *sk = mptcp_to_sock(mptcp);		
		struct tcp_sock *tp = mptcp->tp;
		struct olsched_priv *ol_p = olsched_get_priv(tp);

		all_ended = all_ended && ol_current_send_interval_ended(ol_p->intervals_data[0], tp);
	}

	return all_ended;
}


/* Have we accounted for (acked or lost) enough of the packets that we sent to
 * calculate summary statistics?
 */
bool receive_interval_ended(struct ol_interval_MAB *interval)
{
	if(interval->rcv_ended)
		return true;

	if(interval->snd_ended && !before(interval->known_seq, interval->snd_seq_end)){
		interval->rcv_ended = true;
		return true;
	}

	return false;
}


/* for MAB Solution, check if all intervals in each tp is ended */
bool all_receive_interval_ended(struct sock *meta_sk)
{
	struct tcp_sock *meta_tp = tcp_sk(meta_sk);	
	struct olsched_cb *ol_cb = olsched_get_cb(meta_tp);
	struct ol_interval_MAB *meta_interval = ol_cb->meta_interval[0];
	struct mptcp_tcp_sock *mptcp;
	// struct ol_global *global = olsched_get_cb(meta_tp);
	struct mptcp_cb *mpcb = meta_tp->mpcb;
	bool all_ended = true;

	mptcp_for_each_sub(mpcb, mptcp) {
		struct sock *sk = mptcp_to_sock(mptcp);		
		struct tcp_sock *tp = mptcp->tp;
		struct olsched_priv *ol_p = olsched_get_priv(tp);

		/* TODO valid check, continue */

		all_ended = all_ended && receive_interval_ended(&ol_p->intervals_data[0]);
	}

	return all_ended && receive_interval_ended(meta_interval);;
}

/* Set the red_ratio based on the currently-sending interval
 * and update some interval states.
 */
void start_current_send_interval(struct tcp_sock *tp)
{
	struct olsched_priv *ol_p = olsched_get_priv(tp);
	struct ol_global_MAB *global = ol_p->global_data;
	struct ol_interval_MAB *interval = &ol_p->intervals_data[0];

	interval->snd_time_begin = tp->tcp_mstamp;
	interval->snd_seq_begin = tp->snd_nxt;
	interval->pkts_out_begin = tp->data_segs_out; /* maybe? */
	interval->known_seq = interval->snd_seq_begin; /* init known_seq as the next unacked seq, should be less than snd_next*/
	interval->lost_bytes = 0;
	interval->lost_begin = tp->lost;

	interval->snd_ended = false;
	interval->rcv_ended = false;
	/* when a interval is active, the global ratio is set to the interval's ratio */
	global->red_ratio = interval->red_ratio; 

	printk(KERN_DEBUG "ytxing: tp:%p start interval snd_idx: %u ratio: %u/1024\n", tp, interval->index, global->red_ratio >> 3);
}

void ol_update_arm_probality(struct sock *sk)
{
	struct tcp_sock *tp = tcp_sock(sk);
	struct olsched_priv *ol_p = olsched_get_priv(tp);
	struct ol_interval_MAB *interval = ol_p->intervals_data[0];
	// struct ol_global_MAB *global = ol_p->global_data;

	u32 arm_weight_sum = 0, gamma_over_K = (OLSCHED_UNIT * OLSCHED_GAMMA_MAB) / (OLSCHED_ARMS_NUM * OLSCHED_GAMMA_MAB_BASE);
	u64 temp;
	int i;

	for (i = 0; i < OLSCHED_ARMS_NUM; i++){
		arm_weight_sum += interval->arm_weight[i];
	}

	for (i = 0; i < OLSCHED_ARMS_NUM; i++){
		temp = (interval->arm_weight[i] * OLSCHED_UNIT) / arm_weight_sum;
		temp *= OLSCHED_GAMMA_MAB_BASE - OLSCHED_GAMMA_MAB;
		temp /= OLSCHED_GAMMA_MAB_BASE;
		interval->arm_probability[i] = temp + gamma_over_K;
	}

}

u8 pull_the_arm_accordingly(struct sock *sk)
{
	struct tcp_sock *tp = tcp_sock(sk);
	struct olsched_priv *ol_p = olsched_get_priv(tp);
	struct ol_interval_MAB *interval = ol_p->intervals_data[0];
	// struct ol_global_MAB *global = ol_p->global_data;
	int arm_idx;
	u16 random, probability_cumulated = 0;
	get_random_bytes(&random, OLSCHED_SCALE - 1);
	if (random > OLSCHED_UNIT)
		printk(KERN_DEBUG "ytxing: tp:%p BUG random > OLSCHED_UNIT\n", tp, interval->index, global->red_ratio >> 3);

	for (arm_idx = 0; arm_idx < OLSCHED_ARMS_NUM; arm_idx++){
		probability_cumulated += interval->arm_probability[arm_idx];
		if (random <= probability_cumulated){
			break;
		}
	}

	return arm_idx;
}


/* The Exp3 Algorithm: 
 * 1. set arm probability according to arm weight
 * 2. draw arm randomly according to the probabilities (curr_arm_idx)
 * 3. receive reward of the previous pulled arm (ol_exp3_decide start from here)
 * 4. for each arm, update their weights
 *
 *           1->2->3->4->1->2->3->4->1->2...
 * initialization | in this function
 */
static void ol_exp3_update(struct sock *meta_sk, struct sock *sk)
{
	struct tcp_sock *meta_tp = tcp_sk(meta_sk);
	struct tcp_sock *tp = tcp_sock(sk);
	struct olsched_priv *ol_p = olsched_get_priv(tp);
	struct mptcp_cb *mpcb = meta_tp->mpcb;
	struct ol_interval_MAB *interval = ol_p->intervals_data[0];
	struct ol_global_MAB *global = ol_p->global_data;
	u32 reward, new_ratio;
	int i, j;

	/* update arm weight and probability */
	interval->previous_arm_idx = interval->current_arm_idx;

	/* receive reward from the network */
	reward = ol_calc_reward(interval, sk); // reward << OLSCHED_SCALE, actually. Since 0 < reward <= 1. 
	reward *= OLSCHED_UNIT; /* arm_probability is in OLSCHED_UNIT, so expand the dividend */
	do_div(reward, interval->arm_probability[interval->previous_arm_idx]);
	// interval->arm_weight[interval->previous_arm_idx] update;
	reward *= OLSCHED_UNIT * OLSCHED_GAMMA_MAB;
	do_div(reward, OLSCHED_ARMS_NUM * OLSCHED_GAMMA_MAB_BASE);

	interval->arm_weight[interval->previous_arm_idx] *= ol_exp(reward); /* reward << OLSCHED_SCALE, actually */
	/* Since arm_weight is stored in OLSCHED_UNIT, *= will multiply it to OLSCHED_UNIT^2 */
	interval->arm_weight[interval->previous_arm_idx] >>= OLSCHED_SCALE;

	/* update the probility of previous arm accordingly */
	ol_update_arm_probality(sk);

	/* maybe don't pull the arm immediately, 
	 * since only one tp is allowed to change its behavior at the same time.
	 * TODO
	 */
	interval->current_arm_idx = pull_the_arm_accordingly(sk);
	
	start_current_send_interval(tp);
}



static void update_interval_info_snd(struct ol_interval_MAB *interval, struct tcp_sock *tp)
{
	interval->snd_seq_end = tp->snd_nxt;

	if(interval->snd_time_begin == 0)
		interval->snd_time_begin = tp->tcp_mstamp;

	// printk(KERN_DEBUG "ytxing: tp:%p idx: %u snd_seq: %u->%u(%u)\n", tp, interval->index, interval->snd_seq_begin, interval->snd_seq_end, interval->snd_seq_end - interval->snd_seq_begin);
}

static void update_all_subflow_interval_info_snd(struct tcp_sock *meta_tp)
{
	struct mptcp_tcp_sock *mptcp;
	struct mptcp_cb *mpcb = meta_tp->mpcb;

	mptcp_for_each_sub(mpcb, mptcp) {
		struct sock *sk = mptcp_to_sock(mptcp);		
		struct tcp_sock *tp = mptcp->tp;
		struct olsched_priv *ol_p = olsched_get_priv(tp);

		update_interval_info_snd(ol_p->intervals_data[0], tp);
	}

}

/* ytxing:
 * simply update the known_seq of the interval
 * 
 */
static void update_interval_info_rcv(struct ol_interval_MAB *interval, struct tcp_sock *tp)
{
	struct olsched_priv *ol_p = olsched_get_priv(tp);
	
	u32 current_interval_known_seq = interval->known_seq;
	interval->known_seq = tp->snd_una;
	
	/* ytxing: if this interval just start receiving pkts for the first time */
	if (!after(current_interval_known_seq, interval->snd_seq_begin) && after(interval->known_seq, interval->snd_seq_begin)) {
		interval->delivered_begin = tp->delivered;
		interval->lost_begin = tp->lost;
		interval->rcv_time_begin = tp->tcp_mstamp;
	}
	/* ytxing: if the interval has seen all sent packets and finished sending */
	if (interval->snd_ended && before(current_interval_known_seq, interval->snd_seq_end) && !before(interval->known_seq, interval->snd_seq_end)) {
		interval->delivered_end = tp->delivered;
		interval->lost_end = tp->lost;
		interval->rcv_time_end = tp->tcp_mstamp;
		interval->rcv_ended = true;
		printk(KERN_DEBUG "ytxing: tp:%p (meta:%d) idx:%u RCV END delivered:%u lost:%u known:%u\n", tp, is_meta_tp(tp), interval->index, interval->delivered_end - interval->delivered_end, interval->lost_end - interval->lost_begin, interval->known_seq);
		printk(KERN_DEBUG "ytxing: tp:%p (meta:%d) idx:%u RCV END duration:%llu bandwidth:%llu\n", tp, is_meta_tp(tp), interval->index, interval->rcv_time_end - interval->rcv_time_begin, olsched_get_bandwidth_interval(interval));
		// interval->delivered_end = tp->delivered; /* sure? delivered is not reliable with sack since it might contain some pkts of other intervals*/
	}

}


/* ytxing:	
 * Calculate the redundant quota of this subflow,
 * e.g., how many redundant should it send.
 * force means we go to another interval and refresh the quota accordingly.
 */
static void olsched_check_quota(struct tcp_sock *tp,
					  struct olsched_priv *ol_p, bool force)
{
	struct ol_global_MAB *global = ol_p->global_data;
	// struct ol_interval *curr_interval = NULL;
	// if (global->snd_idx < OLSCHED_INTERVALS_NUM)
	// 	curr_interval = ol_p->intervals_data[global->snd_idx];
	// still has quota, no need to calculate again
	if (!force && (global->red_quota > 0 || global->new_quota > 0))
		return;
	global->red_quota =  (tp->snd_cwnd * global->red_ratio) >> OLSCHED_SCALE;
	global->new_quota = tp->snd_cwnd - global->red_quota;
	printk(KERN_INFO "ytxing: tp:%p use idx:%d red_quota:%u new_quota:%u ratio:%u/1024\n", tp, global->snd_idx, global->red_quota, global->new_quota, global->red_ratio >> 3);
}

/* ytxing:
 * calculate the reward of current interval.
 * When this is called, all tps finish receiving data and get its throughput estimation.
 * reward_j(t) = throughput_meta / sum(throughput_subflow)
 * Then 
 */
static u64 ol_calc_reward(struct ol_interval_MAB *curr_interval, struct sock *curr_sk) {
	struct sock *meta_sk = mptcp_meta_sk(curr_sk);
	struct tcp_sock *meta_tp = tcp_sk(meta_sk), *curr_tp = tcp_sk(curr_sk);	
	struct mptcp_tcp_sock *mptcp;
	struct mptcp_cb *mpcb = meta_tp->mpcb;
	struct olsched_priv *curr_ol_p = olsched_get_priv(curr_tp);
	struct ol_interval_MAB *meta_interval = curr_ol_p->global_data->intervals_data[0];
	struct ol_interval_MAB *interval;
	// struct olsched_cb *ol_cb = olsched_get_cb(meta_tp);

	u64 subflow_throughput_sum = 0, meta_throughput;
	u64 reward;


	meta_throughput = olsched_get_bandwidth_interval(meta_interval);


	/* ytxing:	in this loop, we 
	 *			1. calculate the throughput of all subflows 
	 */
	mptcp_for_each_sub(mpcb, mptcp) {
		struct sock *sk = mptcp_to_sock(mptcp);		
		struct tcp_sock *tp = mptcp->tp;
		struct olsched_priv *ol_p;
		u64 bw_t;

		if (!subflow_is_active(tp) || mptcp_is_def_unavailable(sk)) /* we ignore unavailable subflows*/
			continue;

		ol_p = olsched_get_priv(tp);
		interval = ol_p->intervals_data[0];
		subflow_throughput_sum += olsched_get_bandwidth_interval(interval);
	}
	if (subflow_throughput_sum == 0)
		return 0;
	reward = meta_throughput * OLSCHED_UNIT;
	do_div(reward, subflow_throughput_sum);
	return reward;
}


/* was the ol struct fully inited */
bool ol_valid(struct olsched_priv *ol_p)
{	
	return (ol_p && ol_p->global_data && ol_p->intervals_data && ol_p->global_data->init == 1 && ol_p->intervals_data[0].init == 1);
}

/* Updates the OL model of one sk */
static void ol_process(struct sock *meta_sk, struct sock *sk)
{

	struct tcp_sock *meta_tp = tcp_sk(meta_sk);
	struct tcp_sock *tp = tcp_sk(sk);
	// struct olsched_cb *ol_cb = olsched_get_cb(meta_tp);

	struct olsched_cb *ol_cb = olsched_get_cb(meta_tp);
	struct ol_interval_MAB *meta_interval = ol_cb->meta_interval[0];
	struct olsched_priv *ol_p = olsched_get_priv(tp);
	struct ol_global_MAB *global = ol_p->global_data;
	struct ol_interval_MAB *interval = &ol_p->intervals_data[0]; // ytxing: snd_idx should always be 0 with MAB

	/* maybe some checking? */
	if (!ol_valid(ol_p))
		return;

	/* TODO: maybe return when other tps are processing? */

	/* handle interval in sending state */
	if (global->waiting == false) { 
		// printk(KERN_DEBUG "ytxing: tp:%p update_interval_info_snd\n", tp);
		update_interval_info_snd(interval, tp);
		// update_interval_info_snd(meta_interval, meta_tp);
		if (ol_all_subflow_send_interval_ended(meta_tp)){
			if (!ol_current_send_interval_ended(meta_interval, meta_tp)){
				printk(KERN_DEBUG "ytxing: meta_tp:%p BUG subflows end but not meta\n", meta_tp);
				/* end the sending state this interval */
				meta_interval->pkts_out_end = meta_tp->data_segs_out; /* sure? maybe for non-sack */
				meta_interval->snd_seq_end = meta_tp->snd_nxt;
				meta_interval->snd_time_end = meta_tp->tcp_mstamp;
				meta_interval->snd_ended = true;
			}
			global->waiting = true;
		}

	}

	/* handle interval in receiving state */ 
	update_interval_info_rcv(interval, tp); // TODO update all subflow info
	update_interval_info_rcv(meta_interval, meta_tp);
	// printk(KERN_DEBUG "ytxing: tp:%p update_interval_info_rcv\n", tp);
	/* check if all intervals finish receiving */
	if (all_receive_interval_ended(meta_sk)){
		// printk(KERN_DEBUG "ytxing: tp:%p all_receive_interval_ended\n", tp);
		printk(KERN_DEBUG "ytxing: tp:%p decide\n", tp);
		/* decide new suggested ratio and OL_GLOBAL_STATE*/
		ol_exp3_update(meta_sk, sk); // ytxing: TODO

	}
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
	struct mptcp_tcp_sock *mptcp;
	int active_valid_sks = -1;
	struct sk_buff *skb;
	struct sock *sk, *unavailable_sk, *gambler_sk;
	// int found = 0;
	int i = 0, arm_idx;

	
	/* process ol model immediately, since sending nothing doesn't mean receiving nothing */
	mptcp_for_each_sub(mpcb, mptcp) {
		// printk(KERN_DEBUG "ytxing: tp:%p PROCESS\n", mptcp->tp);
		ol_process(meta_sk, mptcp_to_sock(mptcp));
	}
	gambler_sk = ol_who_is_the_palyer(meta_sk);
	arm_idx = pull_the_arm_accordingly(gambler_sk);

	/* TODO if sk = gambler_sk, set interval according to arm_idx, changing 
	 * if sk != gambler_sk, set the same red ratio, don't move
	 */
	

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
		// printk(KERN_DEBUG "ytxing: i: %d/%d\n", i+1, active_valid_sks);
		sk = get_shorter_rtt_subflow_with_selector(mpcb, unavailable_sk ,&subflow_is_active);
		if (!sk)
			break;
		tp = tcp_sk(sk);
		// printk(KERN_DEBUG "ytxing: better tp:%p srtt: %u\n", tp, (tp->srtt_us >> 3));
		ol_p = olsched_get_priv(tp);
		olsched_correct_skb_pointers(meta_sk, ol_p);

		/* process ol model here? */
		olsched_check_quota(tp, ol_p, false); // ytxing
		if (ol_p->global_data->red_quota >= 1){ // Executing redundant routines
			skb = olsched_next_skb_from_queue(&meta_sk->sk_write_queue,
								ol_p->skb, meta_sk);
		
			if (skb && olsched_use_subflow(meta_sk, active_valid_sks, tp,
							skb)) {
				ol_p->skb = skb;
				ol_p->skb_end_seq = TCP_SKB_CB(skb)->end_seq;
				olsched_update_next_subflow(tp, ol_cb);
				*subsk = (struct sock *)tp;
				// printk(KERN_DEBUG "ytxing skb_len %u (mss?)\n", skb->len);
				ol_p->global_data->red_quota -= 1;
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
				ol_p->global_data->new_quota -= 1;
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
	// printk(KERN_DEBUG "ytxing: tp:%p Nothing to send\n",  tp);

	/* Nothing to send */
	return NULL;
}

static void olsched_release(struct sock *sk)
{
	struct tcp_sock *tp = tcp_sk(sk);
	
	struct olsched_priv *ol_p = olsched_get_priv(tcp_sk(sk));
	struct olsched_cb *ol_cb = olsched_get_cb(tp);

	/* Check if the next subflow would be the released one. If yes correct
	 * the pointer
	 */
	if (ol_cb->next_subflow == tp)
		olsched_update_next_subflow(tp, ol_cb);

	kfree(ol_p->intervals_data);
	kfree(ol_p->global_data);
}

static void olsched_init(struct sock *sk)
{
	struct olsched_priv *ol_p = olsched_get_priv(tcp_sk(sk));
	int i;
	struct olsched_cb *ol_cb = olsched_get_cb(tcp_sk(mptcp_meta_sk(sk)));

	ol_cb->meta_interval = kzalloc(sizeof(struct ol_interval_MAB) * OLSCHED_INTERVALS_NUM * 2,
				 GFP_KERNEL);

	ol_p->intervals_data = kzalloc(sizeof(struct ol_interval_MAB) * OLSCHED_INTERVALS_NUM * 2,
				 GFP_KERNEL);
	ol_p->global_data = kzalloc(sizeof(struct ol_global_MAB),
				 GFP_KERNEL);
	ol_p->global_data->intervals_data = kzalloc(sizeof(struct ol_interval_MAB) * OLSCHED_INTERVALS_NUM * 2,
				 GFP_KERNEL);
	ol_p->global_data->red_ratio = OLSCHED_INIT_RED_RATIO;
	ol_p->global_data->state = OL_START;
	ol_p->global_data->waiting = false;
	ol_p->global_data->red_quota = 0;
	ol_p->global_data->new_quota = 0;

	for (i = 0; i < OLSCHED_INTERVALS_NUM; i++){
		ol_p->intervals_data[i].index = i;
		ol_p->intervals_data[i].init = 1;
		ol_p->intervals_data[i].interval_duration = OLSCHED_INTERVALS_MIN_DURATION;
	}

	for (i = 0; i < OLSCHED_ARMS_NUM; i++){
		ol_p->intervals_data[0]->arm_weight[i] = 1;
	}
	ol_update_arm_probality(sk);

	/* should I pull an arm now? */
	
	ol_p->global_data->init = 1;


	printk(KERN_DEBUG "ytxing: olsched_init tp:%p\n", tcp_sk(sk));
	if(!tcp_is_sack(tcp_sk(sk)))
		printk(KERN_DEBUG "ytxing: olsched_init tp:%p is not sack, maybe bad.\n", tcp_sk(sk));
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
