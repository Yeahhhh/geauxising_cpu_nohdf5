#include <stdlib.h>
#include <math.h>
#include "sim.h"

void
init_temp (PROB_DATATYPE prob[NBETA_PER_WORD][NPROB_MAX])
{
  for (int i = 0; i < NPROB_MAX; ++i) {
    for (int b = 0; b < NBETA_PER_WORD; b++) {
      //prob[b][i] = 2.0f;
      prob[b][i] = UINT32_MAX;
    }
  }
}



// pre-compute propabilities
// load beta from global memory, compute, save propabilities in shared memory

/*
  bidx        0     1     2     3     4     5     6     7
  energy     -6+H  -6-H  -4+H  -4-H  -2+H  -2-H   0+H   0-H


  even
  energy = bidx - 6 + H
  ---------------------------------------
  bidx        0     2     4     6
  energy     -6+H  -4+H  -2+H   0+H


  odd
  energy = bidx - 7 - H = bidx - 6 + H - (1 + 2H)
  ----------------------------------------------
  bidx        1     3     5     7
  energy     -6-H  -4-H  -2-H   0-H
*/

void
compute_temp
(PROB_DATATYPE prob[NBETA_PER_WORD][NPROB_MAX], float *temp_beta_shared, int word)
{
  // for -2 < H < 2, it is OK to compute only first 8 elements of prob[14]
  // keep the rest unchanged

  for (int i = 0; i < NPROB; ++i) {
    for (int b = 0; b < NBETA_PER_WORD; b++) {
      float mybeta = temp_beta_shared[NBETA_PER_WORD * word + b];
      float energy = i - 6 - H - ((-1 * H * 2.0f + 1.0f) * (i & 1));
      //prob[b][i] = expf (2 * energy * mybeta);
      prob[b][i] = expf (2 * energy * mybeta) * UINT32_MAX;
    }
  }

}





// propose a temperature shuffle (inside a lattice)
// mode 0: 0 swap 1 , 2 swap 3 , 4 swap 5 , ...
// mode 1: 0 , 1 swap 2 , 3 swap 4 , ...
void
shuffle (int *temp_idx, float *temp_beta, int *E, int mode)
{
  int order[NBETA];

  for (int i = 0; i < NBETA; ++i)
    order[temp_idx[i]] = i;

  const int bound = NBETA / 2 - mode;
  for (int i = 0; i < bound; ++i) {
    const int idx0 = order[(i << 1) + mode];
    const int idx1 = order[(i << 1) + 1 + mode];
    const float delta_E = E[idx0] - E[idx1];
    const float delta_beta = temp_beta[idx0] - temp_beta[idx1];
    const float val = expf (delta_E * delta_beta);

    if (myrand () / UINT32_MAX < val) {
      const int tmp0 = temp_idx[idx0];
      temp_idx[idx0] = temp_idx[idx1];
      temp_idx[idx1] = tmp0;

      const float tmp1 = temp_beta[idx0];
      temp_beta[idx0] = temp_beta[idx1];
      temp_beta[idx1] = tmp1;
    }
  }

}



PROB_DATATYPE
myrand ()
{
  //return 666;
  return rand ();
}
