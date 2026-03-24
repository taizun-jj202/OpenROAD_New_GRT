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

// Runtime knobs used by NEWGR to generate diverse route portfolios.
int NEWGR_CAP_MODEL
    = 0;  // 0=tuned default, 1=legacy squeeze, 2=direct-WL, 3=ultra-direct.
float NEWGR_CAP_SCALE = 1.0f;
float NEWGR_RUDY_WEIGHT_SCALE = 1.0f;
float NEWGR_PIN_DENSITY_SCALE = 1.0f;

#define M2_ADJ_MIN 0.68
#define M2_ADJ_MAX 0.98
#define M2_ADJ_MID 4.2
#define M2_ADJ_K 1.6

#define M3_ADJ_MIN 0.70
#define M3_ADJ_MAX 0.98
#define M3_ADJ_MID 4.2
#define M3_ADJ_K 1.6

#define MID_ADJ_MIN 0.92
#define MID_ADJ_MAX 1.00
#define MID_ADJ_MID 5.5
#define MID_ADJ_K 1.2

#define HIGH_ADJ 0.92

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

		if (NEWGR_CAP_MODEL == 1) {
			// Legacy SPRoute-style squeeze profile: prioritize routability.
			if (layerID <= 0)
				adj = 0.90f;
			else if (layerID == 1)
				adj = 0.45f + (0.90f - 0.45f) / (1.0f + exp(2.0f * (rudy - 3.5f)));
			else if (layerID == 2)
				adj = 0.45f + (0.90f - 0.45f) / (1.0f + exp(2.0f * (rudy - 3.5f)));
			else if (layerID >= 3 && layerID <= 4)
				adj = 0.90f;
			else if (layerID >= 5)
				adj = 0.70f;
			if (rudy > 7.0f)
				adj *= 0.90f;
		} else if (NEWGR_CAP_MODEL == 2) {
			// Direct-WL profile: keep more near-hard capacity to shorten trunks.
			if (layerID <= 0)
				adj = 0.97f;
			else if (layerID == 1)
				adj = 0.80f + (1.00f - 0.80f) / (1.0f + exp(1.2f * (rudy - 5.0f)));
			else if (layerID == 2)
				adj = 0.82f + (1.00f - 0.82f) / (1.0f + exp(1.2f * (rudy - 5.0f)));
			else if (layerID >= 3 && layerID <= 4)
				adj = 0.95f + (1.00f - 0.95f) / (1.0f + exp(1.0f * (rudy - 6.0f)));
			else if (layerID >= 5)
				adj = 0.97f;
			if (rudy < 2.5f)
				adj = (adj + 0.12f > 1.0f) ? 1.0f : (adj + 0.12f);
			else if (rudy > 9.0f)
				adj *= 0.95f;
		} else if (NEWGR_CAP_MODEL == 3) {
			// Ultra-direct profile: preserve near-hard capacity except in severe hotspots.
			if (layerID <= 0)
				adj = 1.00f;
			else if (layerID == 1)
				adj = 0.90f + (1.00f - 0.90f) / (1.0f + exp(0.9f * (rudy - 8.5f)));
			else if (layerID == 2)
				adj = 0.92f + (1.00f - 0.92f) / (1.0f + exp(0.9f * (rudy - 8.5f)));
			else if (layerID >= 3 && layerID <= 4)
				adj = 0.95f + (1.00f - 0.95f) / (1.0f + exp(0.8f * (rudy - 9.5f)));
			else if (layerID >= 5)
				adj = 0.98f;

			if (rudy < 6.0f)
				adj = 1.0f;
			else if (rudy > 12.0f)
				adj *= 0.92f;
		} else {
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

			// Keep near-hard capacity in sparse regions to avoid detours.
			if (rudy < 1.5f) {
				adj = (adj + 0.08f > 1.0f) ? 1.0f : (adj + 0.08f);
			} else if (rudy > 7.0f) {
				adj *= 0.90f;
			}
		}

		adj *= NEWGR_CAP_SCALE;
		if (adj < 0.05f)
			adj = 0.05f;
		if (adj > 1.0f)
			adj = 1.0f;

		int adjusted = static_cast<int>((float) x * adj);
		if (adjusted < 1)
			adjusted = 1;
		if (adjusted > x)
			adjusted = x;
		return adjusted;
	}
}

#define MAX_FLOW 0

//#define OPENDB 

#endif
