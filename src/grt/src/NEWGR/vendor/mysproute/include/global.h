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
// Runtime-switchable profile used by NEWGR candidate evaluation.
#define NEWGR_CAP_PROFILE_BALANCED 0
#define NEWGR_CAP_PROFILE_WL_FOCUSED 1
#define NEWGR_CAP_PROFILE_DR_FOCUSED 2
#define NEWGR_CAP_PROFILE_3D_SHORT 3
int newgr_capacity_profile = NEWGR_CAP_PROFILE_BALANCED;

// NEWGR hybrid policy:
// keep SPRoute-style congestion awareness, but bias toward FastRoute-like
// shorter paths by reserving less capacity except in strong hotspots.
#define M2_ADJ_MIN 0.78
#define M2_ADJ_MAX 1.00
#define M2_ADJ_MID 5.0
#define M2_ADJ_K 1.4

#define M3_ADJ_MIN 0.80
#define M3_ADJ_MAX 1.00
#define M3_ADJ_MID 5.0
#define M3_ADJ_K 1.4

#define MID_ADJ_MIN 0.90
#define MID_ADJ_MAX 1.00
#define MID_ADJ_MID 6.0
#define MID_ADJ_K 1.2

#define HIGH_ADJ 0.98

#define OBS_NO_STOP 0 // 1 == go through OBS, 0 == hard stop

#define DEFAULT_RUDY_WEIGHT 0.04
#define NO_RUDY 0
#define DEFAULT_PIN_DENSITY_WEIGHT 1.0

#define SMALL_NET_THRSHD 30


int GLOBAL_CAP_ADJ(int x, float rudy, int layerID) //layerID starting from 0, i.e. 0 = metal1, 1 = metal2, 2 = metal3
{
	if(x == 0)
		return 0;
	else {
		float adj;
		const bool wl_profile = (newgr_capacity_profile == NEWGR_CAP_PROFILE_WL_FOCUSED);
		const bool dr_profile = (newgr_capacity_profile == NEWGR_CAP_PROFILE_DR_FOCUSED);
		const bool short3d_profile
			= (newgr_capacity_profile == NEWGR_CAP_PROFILE_3D_SHORT);
		// In low congestion regions, keep hard capacity to avoid unnecessary
		// global-route detours that tend to increase wirelength.
		if (rudy < 1.0) {
			adj = wl_profile ? 1.0f : (dr_profile ? 0.93f : (short3d_profile ? 0.98f : 0.97f));
		} else if (rudy < 2.0) {
			adj = wl_profile ? 0.98f : (dr_profile ? 0.88f : (short3d_profile ? 0.94f : 0.95f));
		} else

		if(layerID == 1) //metal2
			adj = (float) M2_ADJ_MIN + (float)(M2_ADJ_MAX - M2_ADJ_MIN) / (1.0 + exp(M2_ADJ_K * (rudy - M2_ADJ_MID)));
		else if(layerID == 2) //metal3
			adj = (float) M3_ADJ_MIN + (float)(M3_ADJ_MAX - M3_ADJ_MIN) / (1.0 + exp(M3_ADJ_K * (rudy - M3_ADJ_MID)));
		else if(layerID >= 3 && layerID <= 4) //metal4 metal5
			adj = (float) MID_ADJ_MIN + (float)(MID_ADJ_MAX - MID_ADJ_MIN) / (1.0 + exp(MID_ADJ_K * (rudy - MID_ADJ_MID)));
		else if(layerID >= 5) // metal6 and above
			adj = HIGH_ADJ;
		else
			adj = 1.0;

		if (dr_profile) {
			if (layerID <= 2) {
				adj -= 0.04f;
			} else if (layerID >= 5) {
				adj -= 0.02f;
			}
		}
		// CUGR-like 3D preference for NEWGR:
		// keep lower-layer reserve in congestion, but free upper layers so
		// the maze can climb and avoid long 2D detours.
		if (short3d_profile) {
			if (layerID <= 2 && rudy > 2.0f) {
				adj -= 0.05f;
			} else if (layerID <= 4 && rudy > 4.0f) {
				adj -= 0.03f;
			}
			if (layerID >= 5) {
				adj += (rudy > 3.0f) ? 0.03f : 0.02f;
			}
		}

		// Keep a reserve in severe hotspots; reserve depth depends on profile.
		const float hotspot_adj_floor
			= wl_profile ? 0.92f : (dr_profile ? 0.82f : (short3d_profile ? 0.85f : 0.88f));
		if (rudy > 8.0f && adj > hotspot_adj_floor) {
			adj = hotspot_adj_floor;
		}
		if (adj > 1.0f) {
			adj = 1.0f;
		}
		if (adj < 0.65f) {
			adj = 0.65f;
		}

		//adj = 0.9;
		return ((float) x * adj < 2)?  1 : x * adj;
	}
}

#define MAX_FLOW 0

//#define OPENDB 

#endif
