#ifndef GLOBAL_H
#define GLOBAL_H

#include "defDataBase.h"
#include "lefDataBase.h"

parser::defDataBase defDB;   
parser::lefDataBase lefDB;

int numThreads;
int acc_count;
int n_small_undone;
float max_rudy;

#define M2_ADJ_MIN 0.68
#define M2_ADJ_MAX 0.99
#define M2_ADJ_MID 4.2
#define M2_ADJ_K 1.6

#define M3_ADJ_MIN 0.70
#define M3_ADJ_MAX 0.99
#define M3_ADJ_MID 4.2
#define M3_ADJ_K 1.6

#define MID_ADJ_MIN 0.92
#define MID_ADJ_MAX 1.00
#define MID_ADJ_MID 5.5
#define MID_ADJ_K 1.2

#define HIGH_ADJ 0.93

#define OBS_NO_STOP 0 // 1 == go through OBS, 0 == hard stop

#define DEFAULT_RUDY_WEIGHT 0.04
#define NO_RUDY 0
#define DEFAULT_PIN_DENSITY_WEIGHT 0.75

#define SMALL_NET_THRSHD 30


int GLOBAL_CAP_ADJ(int x, float rudy, int layerID) //layerID starting from 0, i.e. 0 = metal1, 1 = metal2, 2 = metal3
{
	if(x == 0)
		return 0;
	else {
		float adj = 1.0f;

		if(layerID <= 0) //metal1
			adj = 0.95f;
		else if(layerID == 1) //metal2
			adj = (float) M2_ADJ_MIN + (float)(M2_ADJ_MAX - M2_ADJ_MIN) / (1.0 + exp(M2_ADJ_K * (rudy - M2_ADJ_MID)));
		else if(layerID == 2) //metal3
			adj = (float) M3_ADJ_MIN + (float)(M3_ADJ_MAX - M3_ADJ_MIN) / (1.0 + exp(M3_ADJ_K * (rudy - M3_ADJ_MID)));
		else if(layerID >= 3 && layerID <= 4) //metal4 metal5
			adj = (float) MID_ADJ_MIN + (float)(MID_ADJ_MAX - MID_ADJ_MIN) / (1.0 + exp(MID_ADJ_K * (rudy - MID_ADJ_MID)));
		else if(layerID >= 5) // metal6 and above
			adj = HIGH_ADJ;

		// Keep near-hard capacity in sparse regions to avoid unnecessary detours.
		if (rudy < 1.8f) {
			adj = (adj + 0.10f > 1.0f) ? 1.0f : (adj + 0.10f);
		} else if (rudy > 7.0f) {
			adj *= 0.90f;
		}

		const int adjusted = static_cast<int>((float) x * adj);
		return (adjusted < 1) ? 1 : adjusted;
	}
}

#define MAX_FLOW 0

//#define OPENDB 

#endif
