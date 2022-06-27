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
#include <linux/sort.h>


/* ytxing: yue */
typedef __u8 u8;
typedef __u16 u16;
typedef __u64 u64;
typedef __s8 s8;
typedef __s32 s32;
typedef __s64 s64;

static bool DEBUG_FIX_ARM __read_mostly = false; 
module_param(DEBUG_FIX_ARM, bool, 0644);	/*模块参数类型，权限*/
MODULE_PARM_DESC(DEBUG_FIX_ARM, "fix the arm");	/*模块参数说明*/

static int DEBUG_FIXED_ARM_IDX __read_mostly = 2; 
module_param(DEBUG_FIXED_ARM_IDX, int, 0644);	/*模块参数类型，权限*/
MODULE_PARM_DESC(DEBUG_FIXED_ARM_IDX, "fix the arm idx");	/*模块参数说明*/

#define DEBUG_USE_DECOUPLED_BWD false
#define DEBUG_USE_MAX_BWD true

static bool DEBUG_USE_NEW_EPOCH __read_mostly = true; 
module_param(DEBUG_USE_NEW_EPOCH, bool, 0644);	/*模块参数类型，权限*/
MODULE_PARM_DESC(DEBUG_USE_NEW_EPOCH, "use reward monitor?");	/*模块参数说明*/

static bool DEBUG_USE_GAMMA_TUNING __read_mostly = true; 
module_param(DEBUG_USE_GAMMA_TUNING, bool, 0644);	/*模块参数类型，权限*/
MODULE_PARM_DESC(DEBUG_USE_GAMMA_TUNING, "use reward monitor?");	/*模块参数说明*/

static bool DEBUG_USE_UCB __read_mostly = true; 
module_param(DEBUG_USE_UCB, bool, 0644);	/*模块参数类型，权限*/
MODULE_PARM_DESC(DEBUG_USE_UCB, "if true, use UCB; else use Exp3");	/*模块参数说明*/

static bool DEBUG_USE_UCB_MON __read_mostly = true; 
module_param(DEBUG_USE_UCB_MON, bool, 0644);	/*模块参数类型，权限*/
MODULE_PARM_DESC(DEBUG_USE_UCB_MON, "use UCB reward monitor?");	/*模块参数说明*/

static bool DEBUG_USE_ACT_RWD __read_mostly = true; 
module_param(DEBUG_USE_ACT_RWD, bool, 0644);	/*模块参数类型，权限*/
MODULE_PARM_DESC(DEBUG_USE_ACT_RWD, "if true, activate the reward");	/*模块参数说明*/

#define DEBUG_USE_FORCE_DECOUPLE false

#define OLSCHED_INTERVALS_NUM 1
#define OLSCHED_INTERVALS_MIN_DURATION 40 * USEC_PER_MSEC /* minimal duration(us) */

/* ytxing:
 * tried 3/2, it seems that longer duration is better.
 */
#define OLSCHED_INTERVALS_DURATION_N_RTT (4/2) /* n * RTT */
#define OLSCHED_INTERVALS_TIMEOUT OLSCHED_INTERVALS_MIN_DURATION * 3
#define OLSCHED_DCP_BWD_TIMEOUT OLSCHED_INTERVALS_MIN_DURATION * 100 /* a timeout for decoupled bandwidth */
#define OLSCHED_MAX_BWD_TIMEOUT OLSCHED_INTERVALS_MIN_DURATION * 100 /* a timeout for decoupled bandwidth */
#define OLSCHED_INTERVAL_MIN_PACKETS 30

#define OLSCHED_SAFE_MAX_WEIGHT 0x0000000fffffffff /* u32 */
#define OLSCHED_GAME_ROUND 4
#define OLSCHED_MIN_QUOTA 12

#define OLSCHED_SCALE 13
#define OLSCHED_UNIT (1 << OLSCHED_SCALE)
#define OLSCHED_MAX_RED_RATIO OLSCHED_UNIT
#define OLSCHED_MIN_RED_RATIO 1
#define OLSCHED_INIT_RED_RATIO (1 * OLSCHED_UNIT)


/* ytxing: for utility function calculation */
static int OLSCHED_LO_GAMMA_MAB __read_mostly = 20; 
module_param(OLSCHED_LO_GAMMA_MAB, int, 0644);	
MODULE_PARM_DESC(OLSCHED_LO_GAMMA_MAB, "OLSCHED_LO_GAMMA_MAB");	

// #define OLSCHED_LO_GAMMA_MAB 30 /* div by OLSCHED_GAMMA_MAB_BASE */ 
static int OLSCHED_HI_GAMMA_MAB __read_mostly = 50; 
module_param(OLSCHED_HI_GAMMA_MAB, int, 0644);
MODULE_PARM_DESC(OLSCHED_HI_GAMMA_MAB, "OLSCHED_HI_GAMMA_MAB");	

static int OLSCHED_ACTIVATION_B __read_mostly = 4; 
module_param(OLSCHED_ACTIVATION_B, int, 0644);
MODULE_PARM_DESC(OLSCHED_ACTIVATION_B, "OLSCHED_ACTIVATION_B");	
#define OLSCHED_GAMMA_MAB_BASE 100



// static const u16 ol_arm_to_red_ratio[4] = {
// 	OLSCHED_UNIT * 1,
// 	OLSCHED_UNIT * 2 / 3,
// 	OLSCHED_UNIT * 1 / 3,
// 	OLSCHED_MIN_RED_RATIO
// };
static const u16 ol_arm_to_red_ratio[3] = {
	OLSCHED_UNIT * 1,
	OLSCHED_UNIT * 1 / 2,
	OLSCHED_MIN_RED_RATIO
};
// static const u16 ol_arm_to_red_ratio[2] = {
// 	OLSCHED_UNIT * 1,
// 	OLSCHED_MIN_RED_RATIO
// };
#define OLSCHED_ARMS_NUM sizeof(ol_arm_to_red_ratio)/sizeof(ol_arm_to_red_ratio[0])

/* new epoch or gamma */
static int OLSCHED_EPOCH_DURATION_PER_ARM __read_mostly = 25; 
module_param(OLSCHED_EPOCH_DURATION_PER_ARM, int, 0644);	/*模块参数类型，权限*/
MODULE_PARM_DESC(OLSCHED_EPOCH_DURATION_PER_ARM, "OLSCHED_EPOCH_DURATION_PER_ARM");	/*模块参数说明*/
#define OLSCHED_EPOCH_DURATION OLSCHED_ARMS_NUM * OLSCHED_EPOCH_DURATION_PER_ARM
#define OLSCHED_HI_GAMMA_DURATION OLSCHED_ARMS_NUM * OLSCHED_EPOCH_DURATION_PER_ARM

static int OLSCHED_MDEV_MUL_UP __read_mostly = 30; 
module_param(OLSCHED_MDEV_MUL_UP, int, 0644);	/*模块参数类型，权限*/
MODULE_PARM_DESC(OLSCHED_MDEV_MUL_UP, "mul");	/*模块参数说明*/
#define OLSCHED_MDEV_MUL OLSCHED_MDEV_MUL_UP/10
static int OLSCHED_CHANGING_THR __read_mostly = 4; 
module_param(OLSCHED_CHANGING_THR, int, 0644);	/*模块参数类型，权限*/
MODULE_PARM_DESC(OLSCHED_CHANGING_THR, "thr");	/*模块参数说明*/
// #define OLSCHED_CHANGING_THR 4

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
	u32 debug1_begin;
	u32 lost_begin; 		/*pkts*/

	/* status for receiving state */
	u64	rcv_time_begin;
	u64	rcv_time_end;
	u32 known_seq;			/* known send sequence (through (s)ack) */
	u32 delivered_end; 		/*pkts*/
	u32 debug1_end;
	u32 lost_end; 			/*pkts*/
	u32 lost_bytes;

	u32	interval_duration;

	u8	snd_end_ready:1,
		snd_ended:1,
		rcv_ended:1,
		init:1,
		unusued:4;
};

enum OL_MON_STATE {
	OL_CHANGE, // this tp is playing MAB game
	OL_STAY, // another is playing game, I have to wait
};

/* it is to store the smoothed_arm_reward, loss rate and RTT fo a subflow in current "stable" condition */
struct ol_monitor {
	u64 smoothed_arm_reward[OLSCHED_ARMS_NUM]; // update like smoothed rtt
	u64 arm_reward_mdev[OLSCHED_ARMS_NUM]; // update like smoothed rtt
	enum OL_MON_STATE state;
	s8 changing_count[OLSCHED_ARMS_NUM]; /* exceeding OLSCHED_CHANGING_THR indicates that the model need a restart */
	u16 epoch_duration; /* current epoch should not end until this becomes zero */
	u8 hi_gamma_duration; /* how many intervals left to change to a smaller gamma? not sure */
	u64 avg_changed_reward[OLSCHED_ARMS_NUM];
	u8 	hi_gamma_flag:1,
		unused:7;
};

struct ol_global {
	u16	red_quota;	// ytxing: how many redundant pkts should this subflow send
	u16	new_quota;	// ytxing: how many new pkts should this subflow send
	u16	red_ratio;	// ytxing: the ratio to calculate current quota, the same as the red_ratio current interval >> OLSCHED_SCALE

	u32 last_time_delivered; /* pkts, for quota calculation */

	u64 decoupled_bandwdith; /* average bandwidth estimated in "redundant intervals" */
	u64	bwd_update_time;

	/* intervals info */
	u8	snd_idx:4, 	// ytxing: the interval index currently sending pkts (< OLSCHED_INTERVALS_NUM)
		rcv_idx:4;	// ytxing: the interval index currently receiving pkts (< OLSCHED_INTERVALS_NUM)

	u8	waiting:1,	// ytxing: ture if all intervals finish sending but not stop receiving
		moving:1,	// ytxing: ture if olsched move toward one direction, stop when utility decreases
		first_move:1,
		init:1,
		unused:4; 
};

struct ol_gambler {
	u8 previous_arm_idx;
	u8 current_arm_idx;
	u64 arm_weight[OLSCHED_ARMS_NUM];
	u16 arm_probability[OLSCHED_ARMS_NUM];
	u32 curr_gamma;

	u32 arm_count[OLSCHED_ARMS_NUM];
	u8	force_decoupled:1,
		unused:7;
	
	/* for UCB */
	u64 arm_avg_reward[OLSCHED_ARMS_NUM];
	u32 arm_count_total;
};

/* Struct to store the data of a single subflow */
struct ol_priv {
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

/* The MPTCP_SCHED_SIZE limitation */
struct ol_priv_out {
	struct ol_priv *ol_p;
};

/* Struct to store the data of the control block */
struct ol_cb {
	/* The next subflow where a skb should be sent or NULL */
	struct tcp_sock *next_subflow;
	struct ol_interval *meta_interval;
	struct ol_gambler *gambler;
	struct ol_monitor *monitor;
	u32 count;
};

struct ol_cb_out {
	struct ol_cb* ol_cb;
};

/* Returns the socket data from a given subflow socket */
static struct ol_priv_out *ol_get_priv_out(struct tcp_sock *tp)
{
	return (struct ol_priv_out *)&tp->mptcp->mptcp_sched[0];
}

/* Returns the socket data from a given subflow socket */
static struct ol_priv *ol_get_priv(struct tcp_sock *tp)
{
	struct ol_priv_out *ol_p_out = ol_get_priv_out(tp);
	return (struct ol_priv *)ol_p_out->ol_p;
}

/* Returns the control block data from a given meta socket */
static struct ol_cb_out *ol_get_cb_out(struct tcp_sock *tp)
{
	return (struct ol_cb_out *)&tp->mpcb->mptcp_sched[0];
}

/* Returns the control block data from a given meta socket */
static struct ol_cb *ol_get_cb(struct tcp_sock *tp)
{
	struct ol_cb_out *ol_cb_out = ol_get_cb_out(tp);
	return (struct ol_cb *)ol_cb_out->ol_cb;
}

/* get x = number * OLSCHED_UNIT, return (e^number)*OLSCHED_UNIT */
static u64 ol_exp(u32 x)
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

/* get x = number * OLSCHED_UNIT, return (e^b(number-1)*OLSCHED_UNIT */
static u64 ol_activate_reward(u64 x)
{

	u64 e_bx;
	u64 e_b;
	u64 b = OLSCHED_ACTIVATION_B;
	if (! DEBUG_USE_ACT_RWD) {
		return x;
	}
	e_bx = ol_exp(b * x);
	e_b = ol_exp(b * OLSCHED_UNIT);
	e_bx *= OLSCHED_UNIT;
	do_div(e_bx, e_b);
	return min_t(u64, e_bx, x);
}
static int ol_get_active_valid_sks_num(struct sock *meta_sk)
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
		
	// mptcp_debug("ytxing: active_valid_sks: %d/%d\n", active_valid_sks, i);
	return active_valid_sks;
}

/* ytxing:	1. bandwidth = cwnd * mss / rtt 
 *			2. bandwidth = delivered_byte / time_interval
 *			in Bps
 */
static u64 ol_get_bandwidth_interval(struct ol_interval *interval, struct tcp_sock *tp)
{
	u64 bandwidth;
	// u32 mss = tp->mss_cache;
	// u64 byte_rcv = (interval->delivered_end - interval->delivered_begin) * mss;
	u64 byte_rcv = interval->snd_seq_end - interval->snd_seq_begin;
	u64 duration = interval->rcv_time_end - interval->snd_time_begin; 
	if (duration <= 0)
		return 0;
	// if (byte_rcv == 0 && is_meta_tp(tp)){
	// 	byte_rcv = interval->snd_seq_end - interval->snd_seq_begin;
	// }
	bandwidth = byte_rcv * USEC_PER_SEC;
	do_div(bandwidth, duration); /* delivery rate, actually */
	// mptcp_debug("ytxing: tp:%p (meta_tp:%d) byte_rcv:%llu duration:%llu bandwidth:%llu\n", tp, is_meta_tp(tp), byte_rcv, duration, bandwidth);
	return bandwidth; /* Bps */
}

/* ytxing:	update max bandwidth, e.g., the bandwidth estimated in "redundant interval" (arm_idx = 0)
 */
static void ol_update_decoupled_bandwidth_interval(struct ol_gambler *gambler, struct ol_interval *interval, struct tcp_sock *tp)
{
	struct ol_priv *ol_p = ol_get_priv(tp);
	u64 curr_bandwidth = ol_get_bandwidth_interval(interval, tp);

	if (gambler->previous_arm_idx != 0){
		if (tp->tcp_mstamp > ol_p->global_data->bwd_update_time + OLSCHED_DCP_BWD_TIMEOUT){
			gambler->force_decoupled = true;
		}
		return;
	}

	ol_p->global_data->decoupled_bandwdith = curr_bandwidth;
	ol_p->global_data->bwd_update_time = tp->tcp_mstamp;
}

/* ytxing:	update max bandwidth, e.g., the bandwidth estimated in "redundant interval" (arm_idx = 0)
 */
static void ol_update_max_bandwidth_interval(struct ol_gambler *gambler, struct ol_interval *interval, struct tcp_sock *tp)
{
	struct ol_priv *ol_p = ol_get_priv(tp);
	u64 curr_bandwidth = ol_get_bandwidth_interval(interval, tp);
	bool nearly_timeout = tp->tcp_mstamp > ol_p->global_data->bwd_update_time + OLSCHED_MAX_BWD_TIMEOUT / 2;
	bool timeout = tp->tcp_mstamp > ol_p->global_data->bwd_update_time + OLSCHED_MAX_BWD_TIMEOUT;

	/* detect higher bandwidth or expire, update */
	if (curr_bandwidth > ol_p->global_data->decoupled_bandwdith){
		ol_p->global_data->decoupled_bandwdith = curr_bandwidth;
		ol_p->global_data->bwd_update_time = tp->tcp_mstamp;
		return;
	}
	if (nearly_timeout && gambler->current_arm_idx == 0){
		ol_p->global_data->decoupled_bandwdith = curr_bandwidth;
		ol_p->global_data->bwd_update_time = tp->tcp_mstamp;
		return;
	}
	gambler->force_decoupled = timeout;

}

/* ytxing:	update max bandwidth, e.g., the bandwidth estimated in "redundant interval" (arm_idx = 0)
 */
static void ol_update_bandwidth_interval(struct ol_gambler *gambler, struct ol_interval *interval, struct tcp_sock *tp)
{
	if (DEBUG_USE_DECOUPLED_BWD){
		ol_update_decoupled_bandwidth_interval(gambler, interval, tp);
	}
	else if (DEBUG_USE_MAX_BWD){
		ol_update_max_bandwidth_interval(gambler, interval, tp);
	}
	else {
		return;
	}
}

/* tweak the suggested ratio and setup four intervals */
static void ol_setup_intervals_MAB(struct sock *sk, int arm_idx)
{	
	struct tcp_sock *tp = tcp_sk(sk);
	struct ol_priv *ol_p = ol_get_priv(tp);
	struct ol_global *global = ol_p->global_data;
	struct ol_interval *interval = &ol_p->intervals_data[0];
	if (!ol_p || !global || !interval){
		return;
	}
	if (arm_idx >=0 && arm_idx < OLSCHED_ARMS_NUM){
		global->red_ratio = ol_arm_to_red_ratio[arm_idx];
		interval->red_ratio = ol_arm_to_red_ratio[arm_idx];
	}
	
	if (interval->delivered_end <= interval->delivered_begin){
		global->last_time_delivered = interval->delivered_end - interval->delivered_begin;
	} else {
		global->last_time_delivered = 0;
	}

	interval->interval_duration = max_t(u32, OLSCHED_INTERVALS_MIN_DURATION, (tp->srtt_us >> 3) * OLSCHED_INTERVALS_DURATION_N_RTT);

	global->snd_idx = 0;
	global->rcv_idx = 0;
	global->waiting = false;
}

/* Have we sent all the data we need to for this interval? Must have at least
 * the minimum number of packets and should have sent 1 RTT worth of data.
 */
bool ol_current_send_interval_end_ready(struct ol_interval *interval, struct tcp_sock *tp)
{
	u32 interval_duration = interval->interval_duration;

	if (interval->snd_time_begin == 0)
		return false;

	/* not enough sending duration */
	if (interval->snd_time_end - interval->snd_time_begin < interval_duration){
		return false;
	}

	if (interval->snd_end_ready == true)
		return true;
		
	// /* not enough packets out and not timeout */
	// if (packets_sent < OLSCHED_INTERVAL_MIN_PACKETS && !timeout)
	// 	return false;

	/* end the sending state this interval */
	interval->snd_end_ready = true;
	return true;
}

/* all send intervals (of all subflows) ended? */
bool ol_all_subflow_send_interval_ended(struct tcp_sock *meta_tp)
{
	struct mptcp_tcp_sock *mptcp;
	// struct ol_global *global = ol_get_cb(meta_tp);
	struct mptcp_cb *mpcb = meta_tp->mpcb;
	bool all_ended = true;

	mptcp_for_each_sub(mpcb, mptcp) {
		struct tcp_sock *tp = mptcp->tp;
		struct ol_priv *ol_p = ol_get_priv(tp);

		all_ended = all_ended && ol_current_send_interval_end_ready(&ol_p->intervals_data[0], tp);
	}

	if (all_ended){
		mptcp_for_each_sub(mpcb, mptcp) {
			struct tcp_sock *tp = mptcp->tp;
			struct ol_priv *ol_p = ol_get_priv(tp);
			struct ol_interval *interval = &ol_p->intervals_data[0];
			// u32 packets_sent = tp->data_segs_out - interval->pkts_out_begin;
			// u32 interval_duration = interval->interval_duration;

			interval->snd_ended = true;
			// mptcp_debug("ytxing: tp:%p (meta:%d) SND END bytes:%u snd_dur:%llu interval_dur:%u\n",\
			tp, is_meta_tp(tp), interval->snd_seq_end - interval->snd_seq_begin, interval->snd_time_end - interval->snd_time_begin, interval->interval_duration);

		}

	}
	
	return all_ended;
}


bool is_receive_interval_ended(struct ol_interval *interval)
{
	// if(interval->rcv_ended)
	// 	return true;

	// if(interval->snd_ended && !before(interval->known_seq, interval->snd_seq_end)){
	// 	interval->rcv_ended = true;
	// 	return true;
	// }

	// return false;
	return interval->rcv_ended;
}

/* for MAB Solution, check if all intervals in each tp is ended */
bool is_all_receive_interval_ended(struct sock *meta_sk)
{
	struct tcp_sock *meta_tp = tcp_sk(meta_sk);	
	struct ol_cb *ol_cb = ol_get_cb(meta_tp);
	struct ol_interval *meta_interval = &ol_cb->meta_interval[0];
	struct mptcp_tcp_sock *mptcp;
	// struct ol_global *global = ol_get_cb(meta_tp);
	struct mptcp_cb *mpcb = meta_tp->mpcb;
	bool all_ended = true;

	mptcp_for_each_sub(mpcb, mptcp) {	
		struct tcp_sock *tp = mptcp->tp;
		struct ol_priv *ol_p = ol_get_priv(tp);

		all_ended = all_ended && is_receive_interval_ended(&ol_p->intervals_data[0]);
	}

	/* if all subflows rcv ended, try to end meta interval */
	if (all_ended)
		all_ended = all_ended && is_receive_interval_ended(meta_interval);

	return all_ended;
}

/* Set the red_ratio based on the currently-sending interval
 * and update some interval states.
 */
void start_current_send_interval(struct tcp_sock *tp)
{
	struct ol_priv *ol_p;
	struct ol_global *global;
	struct ol_interval *interval;
	if (!is_meta_tp(tp)){
		ol_p = ol_get_priv(tp);
		global = ol_p->global_data;
		interval = &ol_p->intervals_data[0];
	} else {
		interval = &ol_get_cb(tp)->meta_interval[0];
	}
	// interval->snd_time_begin = tp->tcp_mstamp;
	interval->snd_time_begin = 0; // mark a new interval
	interval->snd_time_end = 0; 
	interval->rcv_time_begin = 0; // mark a new interval
	interval->rcv_time_end = 0; // mark a new interval
	
	interval->snd_seq_begin = tp->snd_nxt;
	// mptcp_debug("ytxing: tp:%p start interval snd_seq_begin:%u tp->snd_nxt:%u\n", tp, interval->snd_seq_begin, tp->snd_nxt);
	interval->pkts_out_begin = tp->data_segs_out; /* maybe? */
	interval->known_seq = interval->snd_seq_begin; /* init known_seq as the next unacked seq, should be less than snd_next*/
	interval->lost_bytes = 0;
	interval->lost_begin = 0;

	interval->snd_ended = false;
	interval->snd_end_ready = false;
	interval->rcv_ended = false;
	if (is_meta_tp(tp)){
		// mptcp_debug("ytxing: tp:%p (meta_tp) start interval snd_idx: %u\n", tp, interval->index);
		return;
	}
	
	/* when a interval is active, the global ratio is set to the interval's ratio, then check_quota use global infomation */
	global->waiting = false;
	global->red_ratio = interval->red_ratio; 

	// mptcp_debug("ytxing: tp:%p (sub_tp) start interval snd_idx: %u ratio: %u/1024\n", tp, interval->index, global->red_ratio >> 3);
}

void ol_update_arm_probality(struct sock *meta_sk)
{
	struct tcp_sock *meta_tp = tcp_sk(meta_sk);
	struct ol_cb *ol_cb = ol_get_cb(meta_tp);
	struct ol_gambler *gambler = ol_cb->gambler;

	// struct ol_global *global = ol_p->global_data;

	u64 arm_weight_sum = 0, gamma_over_K, gamma;
	u64 temp;
	int i;

	ol_cb->count++;
	
	mptcp_debug("ytxing: tp:%p (meta_tp) ol_update_arm_probality count:%u\n", meta_tp, ol_cb->count);

	if (ol_cb->monitor->hi_gamma_flag && DEBUG_USE_GAMMA_TUNING){
		gamma = OLSCHED_HI_GAMMA_MAB;
		ol_cb->monitor->hi_gamma_duration --;
		if (ol_cb->monitor->hi_gamma_duration == 0){
			ol_cb->monitor->hi_gamma_flag = 0;
		}
	} else {
		gamma = OLSCHED_LO_GAMMA_MAB;
	}

	if (DEBUG_USE_NEW_EPOCH) {
		if (ol_cb->monitor->epoch_duration > 0){
			ol_cb->monitor->epoch_duration --;
		}
	}

	gamma_over_K = (OLSCHED_UNIT * gamma) / (OLSCHED_ARMS_NUM * OLSCHED_GAMMA_MAB_BASE);

	for (i = 0; i < OLSCHED_ARMS_NUM; i++){
		arm_weight_sum += gambler->arm_weight[i];
	}

	for (i = 0; i < OLSCHED_ARMS_NUM; i++){
		temp = (gambler->arm_weight[i] * OLSCHED_UNIT) / arm_weight_sum;
		temp *= OLSCHED_GAMMA_MAB_BASE - gamma;
		temp /= OLSCHED_GAMMA_MAB_BASE;
		gambler->arm_probability[i] = temp + gamma_over_K;
		// mptcp_debug("ytxing: tp:%p (meta_tp) arm_probability[%d]:%u/1024 count:%u\n", meta_tp, i, gambler->arm_probability[i] >> 3, gambler->arm_count[i]);
		mptcp_debug("ytxing: tp:%p (meta_tp) arm_idx:%d temp:%llu gamma_over_K:%llu hi_gamma left:%u epoch_duration:%u\n", meta_tp, i,temp >> 3, gamma_over_K >> 3, ol_cb->monitor->hi_gamma_duration, ol_cb->monitor->epoch_duration);
	}
	

}

u8 pull_the_arm_exp3(struct sock *meta_sk)
{
	struct tcp_sock *meta_tp = tcp_sk(meta_sk);
	struct ol_cb *ol_cb = ol_get_cb(meta_tp);
	struct ol_gambler *gambler = ol_cb->gambler;
	// struct ol_global *global = ol_p->global_data;
	int arm_idx;
	u16 random, probability_cumulated, random_mask;

	if (DEBUG_FIX_ARM){
		gambler->current_arm_idx = DEBUG_FIXED_ARM_IDX;
		return DEBUG_FIXED_ARM_IDX;
	}

	probability_cumulated = 0; 
	random_mask = OLSCHED_UNIT - 1;
	get_random_bytes(&random, sizeof(random));
	// mptcp_debug("ytxing: tp:%p (meta_tp) random:%x random_mask:%x random:%x\n", meta_tp, random, random_mask, random & random_mask);
	random = random & random_mask;
	if (random > OLSCHED_UNIT)
		mptcp_debug("ytxing: tp:%p (meta_tp) BUG random > OLSCHED_UNIT\n", meta_tp);

	for (arm_idx = 0; arm_idx < OLSCHED_ARMS_NUM; arm_idx++){
		probability_cumulated += gambler->arm_probability[arm_idx];
		if (random <= probability_cumulated){
			break;
		}
	}
	gambler->current_arm_idx = arm_idx;
	gambler->arm_count[arm_idx] ++;
	return arm_idx;
}

u8 pull_the_arm_randomly(struct sock *meta_sk)
{
	struct tcp_sock *meta_tp = tcp_sk(meta_sk);
	struct ol_cb *ol_cb = ol_get_cb(meta_tp);
	struct ol_gambler *gambler = ol_cb->gambler;
	// struct ol_global *global = ol_p->global_data;
	int arm_idx, arm_num = OLSCHED_ARMS_NUM;
	u16 random, probability_cumulated, random_mask;

	if (DEBUG_FIX_ARM){
		gambler->current_arm_idx = DEBUG_FIXED_ARM_IDX;
		return DEBUG_FIXED_ARM_IDX;
	}

	probability_cumulated = 0; 
	random_mask = OLSCHED_UNIT - 1;
	get_random_bytes(&random, sizeof(random));
	// mptcp_debug("ytxing: tp:%p (meta_tp) random:%x random_mask:%x random:%x\n", meta_tp, random, random_mask, random & random_mask);
	random = random & random_mask;
	if (random > OLSCHED_UNIT)
		mptcp_debug("ytxing: tp:%p (meta_tp) BUG random > OLSCHED_UNIT\n", meta_tp);

	for (arm_idx = 0; arm_idx < OLSCHED_ARMS_NUM; arm_idx++){
		probability_cumulated += OLSCHED_UNIT / arm_num;
		if (random <= probability_cumulated){
			break;
		}
	}
	gambler->current_arm_idx = arm_idx;
	gambler->arm_count[arm_idx] ++;
	gambler->arm_count_total ++;
	return arm_idx;
}

#define FIXEDPT_BITS 64
#define FIXEDPT_WBITS (FIXEDPT_BITS - OLSCHED_SCALE)
#define OMIT_STDINT
#include "fixedptc.h" 
/* ytxing:
 * calculate the upper confidence bound of specific arm
 * x_bar + sqrt(2ln_n/n_j), where j is the arm_idx
 * return 0 if this arm is never pulled
 */
static u64 ol_get_UCB(struct sock *meta_sk, int arm_idx) {

	struct tcp_sock *meta_tp = tcp_sk(meta_sk);
	struct ol_cb *ol_cb = ol_get_cb(meta_tp);
	struct ol_gambler *gambler = ol_cb->gambler;
	u64 sqrt_item, UCB_result;
	fixedpt fixedpt_ln_n, fixedpt_sqrt_item;
	u64 reward_bar = gambler->arm_avg_reward[arm_idx]; /* in OLSCHED_UINT */
	s64 n = gambler->arm_count_total;
	s64 n_j = gambler->arm_count[arm_idx];

	if (n_j <= 0) // impossible
		return 0;
	fixedpt_ln_n = fixedpt_ln(fixedpt_fromint(n));
	fixedpt_sqrt_item = fixedpt_div(2 * fixedpt_ln_n, fixedpt_fromint(n_j));
	fixedpt_sqrt_item = fixedpt_sqrt(fixedpt_sqrt_item);
	sqrt_item = fixedpt_sqrt_item >> (FIXEDPT_FBITS - OLSCHED_SCALE); /* in OLSCHED_UNIT, i.e., << OLSCHED_SCALE */
	UCB_result = reward_bar + sqrt_item;
	// mptcp_debug("ytxing: tp:%p (meta_tp) arm_idx:%d UCB_result:%llu reward_bar:%llu sqrt_item:%llu\n", meta_tp, arm_idx, UCB_result >> 3, reward_bar >> 3, sqrt_item >> 3);
	
	return UCB_result;
}

u8 pull_the_arm_UCB(struct sock *meta_sk)
{
	struct tcp_sock *meta_tp = tcp_sk(meta_sk);
	struct ol_cb *ol_cb = ol_get_cb(meta_tp);
	struct ol_gambler *gambler = ol_cb->gambler;
	int arm_idx, hi_arm_idx;
	u64 hi_UCB = 0, curr_UCB = 0;

	for (arm_idx = 0; arm_idx < OLSCHED_ARMS_NUM; arm_idx++){
		/* first of all, all arms should be pulled once */
		if (gambler->arm_count[arm_idx] == 0){
			mptcp_debug("ytxing: tp:%p (meta_tp) arm_idx:%d zero count\n", meta_tp, arm_idx);
			hi_arm_idx = arm_idx;
			break;
		}
		curr_UCB = ol_get_UCB(meta_sk, arm_idx);
		
		mptcp_debug("ytxing: tp:%p (meta_tp) arm_idx:%d curr_UCB:%llu avg_reward:%llu count:%u\n", meta_tp, arm_idx, curr_UCB >> 3, gambler->arm_avg_reward[arm_idx] >> 3, gambler->arm_count[arm_idx]);
		if (curr_UCB > hi_UCB){
			hi_UCB = curr_UCB;
			hi_arm_idx = arm_idx;
		}
	}

	gambler->arm_count[hi_arm_idx] ++;
	gambler->arm_count_total ++;
	gambler->current_arm_idx = hi_arm_idx;
	mptcp_debug("ytxing: tp:%p (meta_tp) arm_count_total:%u hi_arm_idx:%u arm_count:%u\n", meta_tp, gambler->arm_count_total, hi_arm_idx, gambler->arm_count[hi_arm_idx]);
	return hi_arm_idx;
}

static void update_interval_info_snd(struct ol_interval *interval, struct tcp_sock *tp)
{
	if (!tp->snd_nxt){
		mptcp_debug("ytxing: tp:%p BUG !meta_tp->snd_nxt\n", tp);
		return;
	}

	if(interval->snd_time_begin == 0 && after(tp->snd_nxt, interval->snd_seq_begin))
		interval->snd_time_begin = tp->tcp_mstamp;

		
	// mptcp_debug("ytxing: tp:%p update_interval_info_snd snd_seq_begin:%u tp->snd_nxt:%u interval->snd_time_begin:%llu\n", tp, interval->snd_seq_begin, tp->snd_nxt, interval->snd_time_begin);
	/* end the sending state this interval */
	interval->pkts_out_end = tp->data_segs_out; /* sure? maybe for non-sack */
	interval->snd_seq_end = tp->snd_nxt;
	interval->snd_time_end = tp->tcp_mstamp;

	// mptcp_debug("ytxing: tp:%p idx: %u snd_seq:%u->%u(%u)\n", tp, interval->index, interval->snd_seq_begin, interval->snd_seq_end, interval->snd_seq_end - interval->snd_seq_begin);
}

/* ytxing:
 * simply update the known_seq of the interval
 * 
 */
static void update_interval_info_rcv(struct ol_interval *interval, struct tcp_sock *tp)
{
	// struct ol_priv *ol_p = ol_get_priv(tp);
	
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
		if (interval->snd_ended){
			interval->rcv_ended = true;
			// mptcp_debug("ytxing: tp:%p (meta:%d) RCV END bytes:%u srtt(us):%u all_dur:%llu bwd:%llu\n", tp, is_meta_tp(tp), interval->snd_seq_end - interval->snd_seq_begin, tp->srtt_us >> 3, interval->rcv_time_end - interval->snd_time_begin, ol_get_bandwidth_interval(interval, tp));
		}

		return;
	}
	/* ytxing: something wrong here */
	if (!is_meta_tp(tp) && interval->snd_ended && !before(interval->known_seq, interval->snd_seq_end)) {
		interval->delivered_end = tp->delivered;
		interval->lost_end = tp->lost;
		interval->rcv_time_end = tp->tcp_mstamp;
		interval->rcv_ended = true;
		// mptcp_debug("ytxing: tp:%p (meta:%d) idx:%u WTF RCV END delivered:%u lost:%u bytes:%u srtt:%u\n", tp, is_meta_tp(tp), interval->index, interval->delivered_end - interval->delivered_begin, interval->lost_end - interval->lost_begin, interval->snd_seq_end - interval->snd_seq_begin, tp->srtt_us >> 3);
		// mptcp_debug("ytxing: tp:%p (meta:%d) idx:%u WTF RCV END all_duration:%llu bandwidth:%llu\n", tp, is_meta_tp(tp), interval->index, interval->rcv_time_end - interval->snd_time_begin, ol_get_bandwidth_interval(interval, tp));
		return;
	}

}

/* ytxing:	
 * Calculate the redundant quota of this subflow,
 * e.g., how many redundant should it send.
 * force means we go to another interval and refresh the quota accordingly.
 */
static void ol_check_quota(struct tcp_sock *tp,
					  struct ol_priv *ol_p, bool force)
{
	struct ol_global *global = ol_p->global_data;
	u32 total_pkts = max_t(u32, OLSCHED_MIN_QUOTA, global->last_time_delivered);
	// struct ol_interval *curr_interval = NULL;
	// if (global->snd_idx < OLSCHED_INTERVALS_NUM)
	// 	curr_interval = ol_p->intervals_data[global->snd_idx];
	// still has quota, no need to calculate again
	if (!force && (global->red_quota > 0 || global->new_quota > 0))
		return;
	global->red_quota =  (total_pkts * global->red_ratio) >> OLSCHED_SCALE;
	global->new_quota = total_pkts - global->red_quota;
	// mptcp_debug("ytxing: tp:%p use idx:%d red_quota:%u new_quota:%u ratio:%u/1024\n", tp, global->snd_idx, global->red_quota, global->new_quota, global->red_ratio >> 3);
}

/* ytxing:
 * calculate the reward of current interval.
 * When this is called, all tps finish receiving data and get its throughput estimation.
 * reward_j(t) = throughput_meta / sum(throughput_subflow)
 * Then 
 */
static u64 ol_calc_reward(struct sock *meta_sk) {
	struct tcp_sock *meta_tp = tcp_sk(meta_sk);
	struct ol_cb *ol_cb = ol_get_cb(meta_tp);
	struct ol_gambler *gambler = ol_cb->gambler;
	struct ol_interval *meta_interval = ol_cb->meta_interval;
	struct ol_interval *interval;
	struct mptcp_tcp_sock *mptcp;
	struct mptcp_cb *mpcb = meta_tp->mpcb;	

	u64 subflow_throughput_sum = 0, meta_throughput;
	u64 reward;


	meta_throughput = ol_get_bandwidth_interval(meta_interval, meta_tp);
	// mptcp_debug("ytxing: tp:%p (meta_tp) calc reward meta_throughput:%llu\n", meta_tp, meta_throughput);


	/* ytxing:	in this loop, we 
	 *			1. calculate the throughput of all subflows 
	 */
	mptcp_for_each_sub(mpcb, mptcp) {
		struct sock *sk = mptcp_to_sock(mptcp);		
		struct tcp_sock *tp = mptcp->tp;
		struct ol_priv *ol_p;

		if (!subflow_is_active(tp) || mptcp_is_def_unavailable(sk)) /* we ignore unavailable subflows*/
			continue;

		ol_p = ol_get_priv(tp);
		interval = &ol_p->intervals_data[0];


		if (DEBUG_USE_DECOUPLED_BWD || DEBUG_USE_MAX_BWD){
			ol_update_bandwidth_interval(gambler, interval, tp);
			subflow_throughput_sum += ol_p->global_data->decoupled_bandwdith;
		}
		else {
			subflow_throughput_sum += ol_get_bandwidth_interval(interval, tp);
		}
	}
	// mptcp_debug("ytxing: tp:%p (meta_tp) calc reward subflow_throughput_sum:%llu\n", meta_tp, subflow_throughput_sum);

	if (subflow_throughput_sum == 0)
		return 0;
	reward = meta_throughput * OLSCHED_UNIT;
	do_div(reward, subflow_throughput_sum);
	// mptcp_debug("ytxing: tp:%p (meta_tp) arm:%d reward:%llu/1024\n", meta_tp, gambler->previous_arm_idx, reward >> 3);
	return reward;
}

/* TODO if the direction of the outsiders is the same as this arm, do not reset */
void ol_check_reward_difference_UCB1(u64 curr_reward, struct sock *meta_sk)
{
	struct tcp_sock *meta_tp = tcp_sk(meta_sk);
	struct ol_cb *ol_cb = ol_get_cb(meta_tp);
	struct ol_gambler *gambler = ol_cb->gambler;
	struct ol_monitor *monitor = ol_cb->monitor;
	s64 tmp_mdev;
	int i;
	bool outsider = false;
	bool is_best = true, is_worst = true;
	int arm_idx = gambler->previous_arm_idx;
	u64 avg_reward = gambler->arm_avg_reward[arm_idx];

	if (monitor->epoch_duration > 0)
		monitor->epoch_duration --;

	for (i = 0; i < OLSCHED_ARMS_NUM; i++){
		is_best = is_best && (avg_reward >= gambler->arm_avg_reward[i]);
		is_worst = is_worst && (avg_reward <= gambler->arm_avg_reward[i]);
	}
	
	tmp_mdev = curr_reward - avg_reward;

	if (monitor->arm_reward_mdev[arm_idx] == 0){
		monitor->arm_reward_mdev[arm_idx] = abs(tmp_mdev);
		mptcp_debug("ytxing: tp:%p (meta_tp) MON idx:%d curr_reward:%llu avg_reward[%d]:%llu %d * arm_reward_mdev:%llu changing_count:%d\n",\
		meta_tp, arm_idx, curr_reward >> 3, arm_idx, avg_reward >> 3, OLSCHED_MDEV_MUL, monitor->arm_reward_mdev[arm_idx] >> 3, monitor->changing_count[arm_idx]);
		return;
	}

	/* check different_reward */
	if ((curr_reward > avg_reward +  monitor->arm_reward_mdev[arm_idx] * OLSCHED_MDEV_MUL) && !is_best){ 
		/* get higher reward, if this arm is not the best, update */
		if (monitor->changing_count[arm_idx] <= 0){
			monitor->changing_count[arm_idx] = 0;
		}
		monitor->changing_count[arm_idx] ++;
		outsider = true;
	}
	else if ((curr_reward < max_t(s64, avg_reward - monitor->arm_reward_mdev[arm_idx] * OLSCHED_MDEV_MUL, 0)) && !is_worst){
		/* get lower reward, if this arm is not the worst, update  */
		if (monitor->changing_count[arm_idx] >= 0){
			monitor->changing_count[arm_idx] = 0;
		}
		monitor->changing_count[arm_idx] --;
		outsider = true;
	}
	else /* in normal range, changing count set to zero*/
	{
		monitor->changing_count[arm_idx] = 0;
	}

	mptcp_debug("ytxing: tp:%p (meta_tp) MON idx:%d curr_reward:%llu avg_reward[%d]:%llu %d * arm_reward_mdev:%llu changing_count:%d epoch:%u\n",\
	meta_tp, arm_idx, curr_reward >> 3, arm_idx, avg_reward >> 3, OLSCHED_MDEV_MUL, monitor->arm_reward_mdev[arm_idx] >> 3, monitor->changing_count[arm_idx], monitor->epoch_duration);

	if (abs(monitor->changing_count[arm_idx]) >= OLSCHED_CHANGING_THR){
		// monitor->state = OL_CHANGE;
		mptcp_debug("ytxing: tp:%p (meta_tp) MON exceed threshold changing_count:%d\n", meta_tp, monitor->changing_count[arm_idx]);

		/* reset the weight of each arm */
		if (DEBUG_USE_UCB_MON){
			mptcp_debug("ytxing: tp:%p (meta_tp) MON new epoch, reset average reward and count\n", meta_tp);
			for (i = 0; i < OLSCHED_ARMS_NUM; i++){
				gambler->arm_count[i] = 0;
				gambler->arm_count_total = 0;
				gambler->arm_avg_reward[i] = 0;
			}
			monitor->changing_count[arm_idx] = 0;
			gambler->arm_count[arm_idx] ++;
			gambler->arm_count_total ++;
			monitor->epoch_duration = OLSCHED_EPOCH_DURATION;
		}
	}


	if (outsider && monitor->epoch_duration == 0){ /* if so, do not update */
		mptcp_debug("ytxing: tp:%p (meta_tp) MON ignore an outsider\n", meta_tp);
		return;
	}

	monitor->arm_reward_mdev[arm_idx] = monitor->arm_reward_mdev[arm_idx] - (monitor->arm_reward_mdev[arm_idx] >> 3) + (abs(tmp_mdev) >> 3);

}

/* The UCB1 Algorithm: 
 * 0. pull each arm once
 * 1. pull the arm i with the highest UCB of each arm (ol_get_UCB)
 * 2. get reward, update the arm_avg_reward of that arm i
 *
 *              0->1->2->1->2->1->2...
 * initialization | in this function
 *
 * this function update the MAB model
 */
static void ol_update_UCB1(struct sock *meta_sk)
{
	struct tcp_sock *meta_tp = tcp_sk(meta_sk);
	struct ol_cb *ol_cb = ol_get_cb(meta_tp);
	struct ol_gambler *gambler = ol_cb->gambler;
	u64 avg_reward, curr_reward;
	gambler->previous_arm_idx = gambler->current_arm_idx;

	// TODO Maybe activate the curr_reward?
	curr_reward = ol_calc_reward(meta_sk);

	ol_check_reward_difference_UCB1(curr_reward, meta_sk);

	/* update the average reward of this arm */
	avg_reward = gambler->arm_avg_reward[gambler->previous_arm_idx];
	avg_reward *= gambler->arm_count[gambler->previous_arm_idx] - 1; // current pulling action is counted, thus minus one.
	avg_reward += curr_reward; // reward << OLSCHED_SCALE, actually. Since 0 < reward <= 1. 
	do_div(avg_reward, gambler->arm_count[gambler->previous_arm_idx]);
	gambler->arm_avg_reward[gambler->previous_arm_idx] = avg_reward;
	
	mptcp_debug("ytxing: tp:%p (meta_tp) idx:%d curr_reward:%llu avg_reward:%llu\n", meta_tp, gambler->previous_arm_idx, curr_reward >> 3, avg_reward >> 3);
}

/* TODO if the direction of the outsiders is the same as this arm, do not reset */
void ol_check_reward_difference(u64 curr_reward, struct sock *meta_sk)
{	
	struct tcp_sock *meta_tp = tcp_sk(meta_sk);
	struct ol_cb *ol_cb = ol_get_cb(meta_tp);
	struct ol_gambler *gambler = ol_cb->gambler;
	struct ol_monitor *monitor = ol_cb->monitor;
	s64 tmp_mdev;
	int i;
	bool outsider = false;
	bool is_best = true, is_worst = true;

	int arm_idx = gambler->previous_arm_idx;

	for (i = 0; i < OLSCHED_ARMS_NUM; i++){
		is_best = is_best && (gambler->arm_weight[arm_idx] >= gambler->arm_weight[i]);
		is_worst = is_worst && (gambler->arm_weight[arm_idx] <= gambler->arm_weight[i]);
	}

	if (curr_reward > OLSCHED_UNIT){
		// mptcp_debug("ytxing: tp:%p (meta_tp) MON curr_reward:%llu too large, arm:%d\n", meta_tp, curr_reward >> 3, arm_idx);
		curr_reward = OLSCHED_UNIT;
	}

	if (monitor->smoothed_arm_reward[arm_idx] == 0){
		monitor->smoothed_arm_reward[arm_idx] = curr_reward;
		monitor->arm_reward_mdev[arm_idx] = 0;	
		mptcp_debug("ytxing: tp:%p (meta_tp) MON curr_reward:%llu smoothed_arm_reward[%d]:%llu %d * arm_reward_mdev:%llu changing_count:%d\n",\
		meta_tp, curr_reward >> 3, arm_idx, monitor->smoothed_arm_reward[arm_idx] >> 3, OLSCHED_MDEV_MUL, monitor->arm_reward_mdev[arm_idx] >> 3, monitor->changing_count[arm_idx]);
		return;
	}
	if (monitor->arm_reward_mdev[arm_idx] == 0){
		tmp_mdev = curr_reward - monitor->smoothed_arm_reward[arm_idx];
		monitor->arm_reward_mdev[arm_idx] = abs(tmp_mdev);
		monitor->smoothed_arm_reward[arm_idx] = monitor->smoothed_arm_reward[arm_idx] - (monitor->smoothed_arm_reward[arm_idx] >> 3) + (curr_reward >> 3);
		mptcp_debug("ytxing: tp:%p (meta_tp) MON curr_reward:%llu smoothed_arm_reward[%d]:%llu %d * arm_reward_mdev:%llu changing_count:%d\n",\
		meta_tp, curr_reward >> 3, arm_idx, monitor->smoothed_arm_reward[arm_idx] >> 3, OLSCHED_MDEV_MUL, monitor->arm_reward_mdev[arm_idx] >> 3, monitor->changing_count[arm_idx]);
		return;
	}

	/* check different_reward */
	if ((curr_reward > monitor->smoothed_arm_reward[arm_idx] +  monitor->arm_reward_mdev[arm_idx] * OLSCHED_MDEV_MUL) && !is_best){ 
		/* get higher reward, if this arm is not the best, update */
		if (monitor->changing_count[arm_idx] <= 0){
			monitor->changing_count[arm_idx] = 0;
		}
		monitor->changing_count[arm_idx] ++;
		outsider = true;
	}
	else if ((curr_reward < max_t(s64, monitor->smoothed_arm_reward[arm_idx] - monitor->arm_reward_mdev[arm_idx] * OLSCHED_MDEV_MUL, 0)) && !is_worst){
		/* get lower reward, if this arm is not the worst, update  */
		if (monitor->changing_count[arm_idx] >= 0){
			monitor->changing_count[arm_idx] = 0;
		}
		monitor->changing_count[arm_idx] --;
		outsider = true;
	}
	else /* in normal range, changing count set to zero*/
	{
		monitor->changing_count[arm_idx] = 0;
	}

	mptcp_debug("ytxing: tp:%p (meta_tp) MON curr_reward:%llu smoothed_arm_reward[%d]:%llu %d * arm_reward_mdev:%llu changing_count:%d\n",\
	meta_tp, curr_reward >> 3, arm_idx, monitor->smoothed_arm_reward[arm_idx] >> 3, OLSCHED_MDEV_MUL, monitor->arm_reward_mdev[arm_idx] >> 3, monitor->changing_count[arm_idx]);

	if (abs(monitor->changing_count[arm_idx]) >= OLSCHED_CHANGING_THR){
		// monitor->state = OL_CHANGE;
		mptcp_debug("ytxing: tp:%p (meta_tp) MON exceed threshold changing_count:%d\n", meta_tp, monitor->changing_count[arm_idx]);

		/* reset the weight of each arm */
		if (DEBUG_USE_NEW_EPOCH && monitor->epoch_duration == 0){
			mptcp_debug("ytxing: tp:%p (meta_tp) MON new epoch, reset weights and counts\n", meta_tp);
			for (i = 0; i < OLSCHED_ARMS_NUM; i++){
				gambler->arm_weight[i] = OLSCHED_UNIT;
				gambler->arm_count[i] = 0;
				monitor->smoothed_arm_reward[i] = 0; /* reset the smoothed_arm_reward */
			}
			monitor->epoch_duration = OLSCHED_EPOCH_DURATION;
		}
		/* tune gamma for a period of time */
		if (DEBUG_USE_GAMMA_TUNING && monitor->hi_gamma_duration == 0){
			mptcp_debug("ytxing: tp:%p (meta_tp) MON hi gamma\n", meta_tp);
			monitor->hi_gamma_flag = 1;
			monitor->hi_gamma_duration = OLSCHED_EPOCH_DURATION; // try
		}

		monitor->changing_count[arm_idx] = 0;
	}

	if (outsider && monitor->epoch_duration == 0){ /* if so, do not update */
		mptcp_debug("ytxing: tp:%p (meta_tp) MON outsider\n", meta_tp);
		return;
	}

	tmp_mdev = curr_reward - monitor->smoothed_arm_reward[arm_idx];
	monitor->arm_reward_mdev[arm_idx] = monitor->arm_reward_mdev[arm_idx] - (monitor->arm_reward_mdev[arm_idx] >> 3) + (abs(tmp_mdev) >> 3);
	monitor->smoothed_arm_reward[arm_idx] = monitor->smoothed_arm_reward[arm_idx] - (monitor->smoothed_arm_reward[arm_idx] >> 3) + (curr_reward >> 3);
	return;
}


/* The Exp3 Algorithm: 
 * 1. set arm probability according to arm weight
 * 2. draw arm randomly according to the probabilities (curr_arm_idx)
 * 3. receive reward of the previous pulled arm (ol_exp3_decide start from here)
 * 4. for each arm, update their weights
 *
 *           1->2->3->4->1->2->3->4->1->2...
 * initialization | in this function
 *
 * this function update the MAB model of the gambler sk
 */
static void ol_update_exp3(struct sock *meta_sk)
{
	struct tcp_sock *meta_tp = tcp_sk(meta_sk);
	struct ol_cb *ol_cb = ol_get_cb(meta_tp);
	struct ol_gambler *gambler = ol_cb->gambler;
	u64 reward, active_reward, exp_reward;
	u64 temp;
	int i;

	/* update arm weight and probability */
	gambler->previous_arm_idx = gambler->current_arm_idx;

	/* receive reward from the network */
	reward = ol_calc_reward(meta_sk); // reward << OLSCHED_SCALE, actually. Since 0 < reward <= 1. 
	active_reward = ol_activate_reward(reward);
	mptcp_debug("ytxing: tp:%p (meta_tp) reward:%llu active_reward:%llu\n", meta_tp, reward >> 3, active_reward >> 3);
	reward = active_reward;
	
	ol_check_reward_difference(reward, meta_sk);
	// reward = monitor->smoothed_arm_reward[gambler->previous_arm_idx];
	// ol_weight_mismatches_reward(meta_sk);

	// check_avg_reward(reward, gambler, monitor);

	reward *= OLSCHED_UNIT; /* arm_probability is in OLSCHED_UNIT, so expand the dividend */
	if (gambler->arm_probability[gambler->previous_arm_idx] == 0){ /* impossible... */
		ol_update_arm_probality(meta_sk);
	}
	do_div(reward, gambler->arm_probability[gambler->previous_arm_idx]);
	// interval->arm_weight[interval->previous_arm_idx] update;
	reward *= OLSCHED_LO_GAMMA_MAB;
	do_div(reward, OLSCHED_ARMS_NUM * OLSCHED_GAMMA_MAB_BASE);
	
	exp_reward = ol_exp(reward);
	// mptcp_debug("ytxing: tp:%p (meta_tp) reward:%llu exp_reward:%llu\n", meta_tp, reward >> 3, exp_reward >> 3);
	// mptcp_debug("ytxing: tp:%p (meta_tp) arm_weight[%d]:%llu\n", meta_tp, gambler->previous_arm_idx, gambler->arm_weight[gambler->previous_arm_idx]);
	temp = gambler->arm_weight[gambler->previous_arm_idx] * exp_reward;
	// mptcp_debug("ytxing: tp:%p (meta_tp) temp:%llu\n", meta_tp, temp);
	temp >>= OLSCHED_SCALE;
	// mptcp_debug("ytxing: tp:%p (meta_tp) temp:%llu\n", meta_tp, temp);

	while (temp > OLSCHED_SAFE_MAX_WEIGHT){	
		
		mptcp_debug("ytxing: tp:%p (meta_tp) shrink all arm_weight\n", meta_tp);
		/********************************************\
		 * in Exp3, arm_weight will exponential grow *
		 * and all arm_weight needs to be shrinked   *
		 * sometime.                                 *
		\********************************************/
		for (i = 0; i < OLSCHED_ARMS_NUM; i++){
			if (!gambler || !gambler->arm_weight)
				continue;
			gambler->arm_weight[i] = max_t(u64, OLSCHED_UNIT, gambler->arm_weight[i] >> 8);
			mptcp_debug("ytxing: tp:%p (meta_tp) shrink arm_weight[%d]:%llu\n", meta_tp, i, gambler->arm_weight[i]);
		}
		temp = (gambler->arm_weight[gambler->previous_arm_idx] * exp_reward) >> OLSCHED_SCALE;
	}

	gambler->arm_weight[gambler->previous_arm_idx] = temp;


	/* Since arm_weight is stored in OLSCHED_UNIT, *= will multiply it to OLSCHED_UNIT^2 */
	// mptcp_debug("ytxing: tp:%p (meta_tp) arm_weight[%d]:%llu\n", meta_tp, gambler->previous_arm_idx, gambler->arm_weight[gambler->previous_arm_idx]);

	/* update the probility of previous arm accordingly */
	ol_update_arm_probality(meta_sk);
}

static void ol_update(struct sock *meta_sk){
	if (DEBUG_USE_UCB){
		ol_update_UCB1(meta_sk);
	}
	else{
		ol_update_exp3(meta_sk);		
	}
}

/* was the ol struct fully inited */
bool ol_valid(struct ol_priv *ol_p)
{	
	return (ol_p && ol_p->global_data && ol_p->intervals_data && ol_p->intervals_data[0].init && ol_p->global_data->init);
}

/* Updates the OL model of one sk */
static void ol_process_all_subflows(struct sock *meta_sk, bool *new_interval_started)
{
	struct tcp_sock *meta_tp = tcp_sk(meta_sk);
	struct sock *sk;
	struct tcp_sock *tp;
	// struct ol_cb *ol_cb = ol_get_cb(meta_tp);

	struct ol_cb *ol_cb = ol_get_cb(meta_tp);
	struct ol_interval *meta_interval = &ol_cb->meta_interval[0];
	struct ol_priv *ol_p;
	struct ol_global *global;
	struct ol_gambler *gambler = ol_cb->gambler;
	struct ol_interval *interval; // ytxing: snd_idx should always be 0 with MAB
	
	struct mptcp_tcp_sock *mptcp;
	struct mptcp_cb *mpcb = meta_tp->mpcb;

	u32 packets_sent;
	u32 interval_duration;
	u8 arm_idx;

	// bool new_gamma_flag;

	*new_interval_started = false;

	mptcp_for_each_sub(mpcb, mptcp) {
		
		sk = mptcp_to_sock(mptcp);
		tp = tcp_sk(sk);
		ol_p = ol_get_priv(tp);
		interval = &ol_p->intervals_data[0];

		if (mptcp_is_def_unavailable(sk)){
			continue;
		}
		if (!ol_valid(ol_p)){
			continue;
		}
		
		/* SND INFO */
		if (meta_interval->snd_ended == false){
			update_interval_info_snd(interval, tp);
			ol_current_send_interval_end_ready(interval, tp);
		}

		/* RCV INFO */
		// if (interval->rcv_ended == false)
		// 	update_interval_info_rcv(interval, tp);
	}

	/* SND INFO META */
	if (meta_interval->snd_ended == false){
		update_interval_info_snd(meta_interval, meta_tp);
		if (ol_all_subflow_send_interval_ended(meta_tp)){
			if (!ol_current_send_interval_end_ready(meta_interval, meta_tp)){
				// mptcp_debug("ytxing: meta_tp:%p BUG subflows snd end but not meta\n", meta_tp);
			}
			/* end the meta snd interval, same as in ol_all_subflow_send_interval_ended */
			packets_sent = meta_tp->data_segs_out - meta_interval->pkts_out_begin;
			interval_duration = meta_interval->interval_duration;

			meta_interval->snd_ended = true;
			// mptcp_debug("ytxing: tp:%p (meta:%d) SND END packets_sent:%u\n", meta_tp, is_meta_tp(meta_tp), packets_sent);
		}
	}
	
	mptcp_for_each_sub(mpcb, mptcp) {
		
		sk = mptcp_to_sock(mptcp);
		tp = tcp_sk(sk);
		ol_p = ol_get_priv(tp);
		interval = &ol_p->intervals_data[0];

		if (mptcp_is_def_unavailable(sk)){
			continue;
		}
		if (!ol_valid(ol_p)){
			continue;
		}
		
		// /* SND INFO */
		// if (meta_interval->snd_ended == false){
		// 	update_interval_info_snd(interval, tp);
		// 	ol_current_send_interval_end_ready(interval, tp);
		// }

		/* RCV INFO */
		if (interval->rcv_ended == false)
			update_interval_info_rcv(interval, tp);
	}

	/* RCV INFO META */
	update_interval_info_rcv(meta_interval, meta_tp);
	if (!is_all_receive_interval_ended(meta_sk)){
		/* if not all intervals finish, return and wait. */
		return;
	}

	if (gambler == NULL)
		return;

	if (ol_get_active_valid_sks_num(meta_sk) >= 2){
		ol_update(meta_sk);
	} else {
		mptcp_debug("ytxing: tp:%p (meta_tp) WTF too few sks, do not update\n", meta_tp);
	}
	
	/* if we need to decouple subflows to estimate bandwidth of them */
	if (DEBUG_USE_FORCE_DECOUPLE && gambler->force_decoupled){
		mptcp_debug("ytxing: tp:%p (meta_tp) force_decoupled\n", meta_tp);
		gambler->force_decoupled = false;
		gambler->current_arm_idx = 0;
		arm_idx = 0;
	} else {
		/* pull an arm for all subflows */
		if (DEBUG_USE_UCB) {
			arm_idx = pull_the_arm_UCB(meta_sk);
		}
		else {
			arm_idx = pull_the_arm_exp3(meta_sk);
		}
		// arm_idx = pull_the_arm_randomly(meta_sk);
	}

	/* setup interval for all subflows */
	mptcp_for_each_sub(mpcb, mptcp) {
		// mptcp_debug("ytxing: tp:%p PROCESS\n", mptcp->tp);
		sk = mptcp_to_sock(mptcp);
		ol_p = ol_get_priv(tcp_sk(sk));
		global = ol_p->global_data;

		if (mptcp_is_def_unavailable(sk)){
			continue;
		}
		if (!ol_valid(ol_p)){
			continue;
		}
		ol_setup_intervals_MAB(sk, arm_idx);

		start_current_send_interval(tcp_sk(sk));
	}
	start_current_send_interval(tcp_sk(meta_sk));
	
	*new_interval_started = true;
}



static bool ol_use_subflow(struct sock *meta_sk,
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
			active_valid_sks = ol_get_active_valid_sks_num(meta_sk);

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

static void ol_update_next_subflow(struct tcp_sock *tp,
					 struct ol_cb *ol_cb)
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
	struct ol_cb *ol_cb = ol_get_cb(meta_tp);
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
			ol_update_next_subflow(tp, ol_cb);
			return (struct sock *)tp;
		}
	}

	mptcp_for_each_sub(mpcb, mptcp) {
		tp = mptcp->tp;

		if (tp == first_tp)
			break;

		if (mptcp_is_available((struct sock *)tp, skb,
				       zero_wnd_test)) {
			ol_update_next_subflow(tp, ol_cb);
			return (struct sock *)tp;
		}
	}

	/* No space */
	return NULL;
}

/* Corrects the stored skb pointers if they are invalid */
static void ol_correct_skb_pointers(struct sock *meta_sk,
					  struct ol_priv *ol_p)
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
static struct sk_buff *ol_next_skb_from_queue(struct sk_buff_head *queue,
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

// ytxing: 	1. find a subflow with minimal RTT;
//			2. get a redundant or normal skb for the subflow
static struct sk_buff *mptcp_ol_next_segment_rtt(struct sock *meta_sk,
					      int *reinject,
					      struct sock **subsk,
					      unsigned int *limit)
{
	struct tcp_sock *meta_tp = tcp_sk(meta_sk);
	struct mptcp_cb *mpcb = meta_tp->mpcb;
	struct ol_cb *ol_cb = ol_get_cb(meta_tp);
	struct ol_priv *ol_p;
	struct tcp_sock *tp;
	int active_valid_sks = -1;
	struct sk_buff *skb;
	struct sock *sk, *unavailable_sk;
	// int found = 0;
	int i = 0;
	bool new_interval_started;

	active_valid_sks = ol_get_active_valid_sks_num(meta_sk);
	
	/* process ol model immediately, since sending nothing doesn't mean receiving nothing */
	if (active_valid_sks >= 2){
		// mptcp_debug("ytxing: ol_process_all_subflows\n");
		ol_process_all_subflows(meta_sk, &new_interval_started);
	}
	

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

	unavailable_sk = NULL;
	/* ytxing: in this loop we find a avaliable sk with rtt as short as possible*/
	for(i = 0; i < active_valid_sks; i++){
		// mptcp_debug("ytxing: i: %d/%d\n", i+1, active_valid_sks);
		sk = get_shorter_rtt_subflow_with_selector(mpcb, unavailable_sk ,&subflow_is_active);
		if (!sk)
			break;
		tp = tcp_sk(sk);
		// mptcp_debug("ytxing: better tp:%p srtt: %u\n", tp, (tp->srtt_us >> 3));
		ol_p = ol_get_priv(tp);
		ol_correct_skb_pointers(meta_sk, ol_p);

		ol_check_quota(tp, ol_p, new_interval_started); // ytxing: if new_interval_started, refresh the quota.
		
		if (ol_p->global_data->red_quota >= 1){ // Executing redundant routines
			skb = ol_next_skb_from_queue(&meta_sk->sk_write_queue,
								ol_p->skb, meta_sk);
		
			if (skb && ol_use_subflow(meta_sk, active_valid_sks, tp,
							skb)) {
				ol_p->skb = skb;
				ol_p->skb_end_seq = TCP_SKB_CB(skb)->end_seq;
				ol_update_next_subflow(tp, ol_cb);
				*subsk = (struct sock *)tp;
				// mptcp_debug("ytxing skb_len %u (mss?)\n", skb->len);
				ol_p->global_data->red_quota -= 1;
				if (TCP_SKB_CB(skb)->path_mask){
					*reinject = -1;
					// mptcp_debug("ytxing: a red pkt\n");
					}
				// mptcp_debug("ytxing: send pkt 1\n");
				return skb;
			}
		}
		else { // send a new skb, no need redundant
			skb = tcp_send_head(meta_sk);
			if (skb && ol_use_subflow(meta_sk, active_valid_sks, tp,
							skb)) {
				ol_p->skb = skb;
				ol_p->skb_end_seq = TCP_SKB_CB(skb)->end_seq;
				ol_update_next_subflow(tp, ol_cb);
				*subsk = (struct sock *)tp;
				// mptcp_debug("ytxing skb_len %u (mss?)\n", skb->len);
				ol_p->global_data->new_quota -= 1;
				if (TCP_SKB_CB(skb)->path_mask)
					*reinject = -1;
				// mptcp_debug("ytxing: send pkt 2\n");	
				return skb;
			}
		}

		// ytxing: 	now this sk is unable to send (maybe cwnd limited)
		//			to a sk with longer rtt.
		unavailable_sk = sk;
	}
	// mptcp_debug("ytxing: tp:%p Nothing to send\n",  tp);

	/* Nothing to send */
	return NULL;
}

static void ol_release(struct sock *sk)
{
	struct tcp_sock *tp = tcp_sk(sk);
	
	struct ol_priv *ol_p = ol_get_priv(tcp_sk(sk));
	struct ol_cb *ol_cb = ol_get_cb(tp);

	/* Check if the next subflow would be the released one. If yes correct
	 * the pointer
	 */
	if (ol_cb->next_subflow == tp)
		ol_update_next_subflow(tp, ol_cb);

	kfree(ol_p->intervals_data);
	kfree(ol_p->global_data);
}

static void ol_init(struct sock *sk)
{
	struct ol_priv *ol_p;
	struct ol_cb *ol_cb;
	int i;
	struct ol_priv_out *ols_p_out = ol_get_priv_out(tcp_sk(sk));
	struct ol_cb_out *ols_cb_out = ol_get_cb_out(tcp_sk(mptcp_meta_sk(sk)));
	ols_p_out->ol_p = kzalloc(sizeof(struct ol_priv), GFP_KERNEL);
	ol_p = ol_get_priv(tcp_sk(sk));
	ols_cb_out->ol_cb = kzalloc(sizeof(struct ol_cb), GFP_KERNEL);
	ol_cb = ol_get_cb(tcp_sk(mptcp_meta_sk(sk)));


	ol_p->intervals_data = kzalloc(sizeof(struct ol_interval) * OLSCHED_INTERVALS_NUM, GFP_KERNEL);
	ol_p->global_data = kzalloc(sizeof(struct ol_global), GFP_KERNEL);
	ol_cb->monitor = kzalloc(sizeof(struct ol_monitor), GFP_KERNEL);

	for (i = 0; i < OLSCHED_INTERVALS_NUM; i++){
		ol_p->intervals_data[i].index = i;
		ol_p->intervals_data[i].init = 1;
		ol_p->intervals_data[i].interval_duration = OLSCHED_INTERVALS_MIN_DURATION;
		ol_p->intervals_data[i].red_ratio = OLSCHED_INIT_RED_RATIO;
		ol_p->intervals_data[i].snd_ended = false;
		ol_p->intervals_data[i].rcv_ended = false;
	}

	ol_p->global_data->red_ratio = OLSCHED_INIT_RED_RATIO;
	ol_p->global_data->last_time_delivered = OLSCHED_MIN_QUOTA;
	ol_p->global_data->waiting = false;
	ol_p->global_data->red_quota = 0;
	ol_p->global_data->new_quota = 0;
	ol_p->global_data->decoupled_bandwdith = 0;
	ol_p->global_data->bwd_update_time = 0;

	ol_cb->count = 0;
	if (!ol_cb->meta_interval){
		ol_cb->meta_interval = kzalloc(sizeof(struct ol_interval) * OLSCHED_INTERVALS_NUM, GFP_KERNEL);
		ol_cb->meta_interval->index = 0;
		ol_cb->meta_interval->init = 1;
		ol_cb->meta_interval->interval_duration = OLSCHED_INTERVALS_MIN_DURATION;
		ol_cb->meta_interval->snd_ended = false;
		ol_cb->meta_interval->rcv_ended = false;
		ol_cb->meta_interval->snd_seq_begin = 0;
		mptcp_debug("ytxing: ol_init meta_interval meta_tp:%p\n", tcp_sk(mptcp_meta_sk(sk)));
	}

	/* if this is the first tp initiated, it will be the first gambler. */
	if (!ol_cb->gambler) {
		ol_cb->gambler = kzalloc(sizeof(struct ol_gambler), GFP_KERNEL);
		
		for (i = 0; i < OLSCHED_ARMS_NUM; i++){
			ol_cb->gambler->arm_weight[i] = OLSCHED_UNIT;
			ol_cb->gambler->arm_count[i] = 0;
			ol_cb->gambler->arm_avg_reward[i] = 0;
		}
		ol_cb->gambler->curr_gamma = OLSCHED_LO_GAMMA_MAB;
		ol_cb->gambler->force_decoupled = false;
		ol_cb->gambler->arm_count_total = 0;
		ol_update_arm_probality(mptcp_meta_sk(sk));
		ol_cb->gambler->current_arm_idx = 0;
		ol_cb->gambler->arm_count[0] ++;
		ol_cb->gambler->arm_count_total ++;
		mptcp_debug("ytxing: ol_init gambler meta_tp:%p\n", tcp_sk(mptcp_meta_sk(sk)));

	}
	
	ol_cb->monitor->state = OL_CHANGE;
	for (i = 0; i < OLSCHED_INTERVALS_NUM; i++){
		ol_cb->monitor->smoothed_arm_reward[i] = 0;
		ol_cb->monitor->arm_reward_mdev[i] = 0;
		ol_cb->monitor->avg_changed_reward[i] = 0;
		ol_cb->monitor->changing_count[i] = 0;
		ol_cb->monitor->epoch_duration = OLSCHED_EPOCH_DURATION;
		ol_cb->monitor->hi_gamma_flag = 1;
		ol_cb->monitor->hi_gamma_duration = OLSCHED_EPOCH_DURATION;
	}

	ol_p->global_data->init = 1;

	start_current_send_interval(tcp_sk(sk));
	if (ol_cb->meta_interval->snd_seq_begin == 0){
		start_current_send_interval(tcp_sk(mptcp_meta_sk(sk))); 
	}


	mptcp_debug("ytxing: ol_init tp:%p\n", tcp_sk(sk));
	mptcp_debug("ytxing: tp:%p DEBUG_USE_NEW_EPOCH:%d DEBUG_USE_GAMMA_TUNING:%d\n", tcp_sk(sk), DEBUG_USE_NEW_EPOCH, DEBUG_USE_GAMMA_TUNING);
	// if(!tcp_is_sack(tcp_sk(sk)))
	// 	mptcp_debug("ytxing: ol_init tp:%p is not sack, maybe bad.\n", tcp_sk(sk));
}

static struct mptcp_sched_ops mptcp_sched_ol = {
	.get_subflow = ol_get_available_subflow,
	.next_segment = mptcp_ol_next_segment_rtt,
	.init = ol_init,
	.release = ol_release,
	.name = "ol",
	.owner = THIS_MODULE,
};

static int __init ol_register(void)
{
	BUILD_BUG_ON(sizeof(struct ol_priv_out) > MPTCP_SCHED_SIZE);
	BUILD_BUG_ON(sizeof(struct ol_cb_out) > MPTCP_SCHED_DATA_SIZE);
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
MODULE_DESCRIPTION("An Online-Learning-Asisted Packet Scheduler(OLAPS) for MPTCP");
MODULE_VERSION("0.21");
