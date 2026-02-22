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

#define M2_ADJ_MIN 0.69
#define M2_ADJ_MAX 1.00
#define M2_ADJ_MID 3.8
#define M2_ADJ_K 2.0

#define M3_ADJ_MIN 0.69
#define M3_ADJ_MAX 1.00
#define M3_ADJ_MID 3.8
#define M3_ADJ_K 2.0

#define MID_ADJ_MIN 0.985
#define MID_ADJ_MAX 1.00
#define MID_ADJ_MID 5.0
#define MID_ADJ_K 2.0

#define HIGH_ADJ 0.93

#define OBS_NO_STOP 0 // 1 == go through OBS, 0 == hard stop

#define DEFAULT_RUDY_WEIGHT 0.1
#define NO_RUDY 0
#define DEFAULT_PIN_DENSITY_WEIGHT 1.0

#define SMALL_NET_THRSHD 30


int GLOBAL_CAP_ADJ(int x, float rudy, int layerID) //layerID starting from 0, i.e. 0 = metal1, 1 = metal2, 2 = metal3
{
	if(x == 0)
		return 0;
	else {
		// Default to no soft-cap reduction for layers that are not explicitly
		// parameterized below (e.g. metal1).
		float adj = 1.0f;

		if(layerID == 1) //metal2
			adj = (float) M2_ADJ_MIN + (float)(M2_ADJ_MAX - M2_ADJ_MIN) / (1.0 + exp(M2_ADJ_K * (rudy - M2_ADJ_MID)));
		else if(layerID == 2) //metal3
			adj = (float) M3_ADJ_MIN + (float)(M3_ADJ_MAX - M3_ADJ_MIN) / (1.0 + exp(M3_ADJ_K * (rudy - M3_ADJ_MID)));
		else if(layerID >= 3 && layerID <= 4) //metal4 metal5
			adj = (float) MID_ADJ_MIN + (float)(MID_ADJ_MAX - MID_ADJ_MIN) / (1.0 + exp(MID_ADJ_K * (rudy - MID_ADJ_MID)));
		else if(layerID >= 5) // metal6 and above
			adj = HIGH_ADJ;
		

		// Keep soft-cap behavior, but avoid overly aggressive truncation on
		// low-capacity edges (e.g. 1.8 -> 1), which can over-constrain routing.
		const float scaled_cap = (float) x * adj;
		const int adjusted_cap = (int) (scaled_cap + 0.5f);
		return (adjusted_cap < 1) ? 1 : adjusted_cap;
	}
}

#define MAX_FLOW 0

//#define OPENDB 

#endif
