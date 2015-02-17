#include "COPYING"


#ifndef SIM_H
#define SIM_H

#include <stdint.h>
#include "parameters.h"

#define PI 3.141592654f

// string length for file names, etc.
#define STR_LENG 64

// Multispin Coding using 32 bit integer
//#include "bit32.h"
#include "bit64.h"

#define NO_OUTPUT



/*
   probablities


   __expf (2 * energy * temp_beta_shared[b]);
   temp_prob_shared[16 * b + (energy << 1) + spin];


   probability = expf (2 * beta * (energy - H * spin))
   prob []

   index energy spin   (energy - H * spin)
   00    -6       -1    -6+H
   01    -6       +1    -6-H
   
   02    -4       -1    -4+H
   03    -4       +1    -4-H
   
   04    -2       -1    -2+H
   05    -2       +1    -2-H
   
   06     0       -1     0+H
   07     0       +1     0-H
   
   08     2       -1     2+H
   09     2       +1     2-H
   
   10     4       -1     4+H
   11     4       +1     4-H
   
   12     6       -1     6+H
   13     6       +1     6-H
 */

#define NPROB 14
#define NPROB_MAX 16
//typedef float PROB_DATATYPE;
typedef uint32_t PROB_DATATYPE;
//#define UINT32_MAX 0xffffffff







// better allign on 4 byte boundary

typedef struct
{
  int E[NBETA][REALIZATION];
  int M[NBETA][REALIZATION];
  float U[NBETA][REALIZATION];		// U = 1 - M4 / (3 * M2 * M2)
} Avrg;


/*
typedef struct
{
  float Q0[NBETA][BD];
  float Qk_real[NBETA][BD];
  float Qk_imag[NBETA][BD];
  float Qk2_real[NBETA][BD];
  float Qk2_imag[NBETA][BD];                                               
} Qk;
*/


typedef struct
{
  //  float e[REALIZATION][NBETA_MAX];
  float q[REALIZATION_HF][NBETA_MAX];
  float qk_real[REALIZATION_HF][NBETA_MAX];
  float qk_imag[REALIZATION_HF][NBETA_MAX];
  float qk2_real[REALIZATION_HF][NBETA_MAX];
  float qk2_imag[REALIZATION_HF][NBETA_MAX];
} St_old;


typedef struct
{
  //  float e[REALIZATION][NBETA_MAX];
  float q[REALIZATION_HF][NBETA_MAX];
  float qk_real[3][REALIZATION_HF][NBETA_MAX];
  float qk_imag[3][REALIZATION_HF][NBETA_MAX];
  float qk2_real[6][REALIZATION_HF][NBETA_MAX];
  float qk2_imag[6][REALIZATION_HF][NBETA_MAX];
} St;


typedef struct
{
  int idx;
  float beta;
} Temp;


typedef struct
{
  //curandState *gpuseed;
  MSC_DATATYPE *lattice;
  MSC_DATATYPE *lattice1;
  Temp *temp;
  St *st;
} Parameter;






// host_func.c
void host_init_J (MSC_DATATYPE * l);
void host_init_S (MSC_DATATYPE * l);
void host_init_lattice (MSC_DATATYPE * l);
void host_save_st (St * st, char *mydir, int node, int device);
void host_save_st_hdf5 (St * st, char *mydir, int node, int device);
void host_init_temp (Temp * temp, float beta_low, float beta_high);
void host_report_speed_title ();
void host_report_speed (double start, double stop, int iter, char *event);
void host_usage (char *bin);
void host_summary (float beta_low, float beta_high, char *mydir);

//host_launcher.c
void host_launcher (float beta_low, float beta_high, char* mydir, int node, int device);

//host_kernel.c
void host_kernel_warmup (MSC_DATATYPE * lattice, Temp * temp);
void host_kernel_swap (MSC_DATATYPE * lattice, Temp * temp, St * st, int rec);




// kernel_l3.c
void init_temp (PROB_DATATYPE prob[NBETA_PER_WORD][NPROB_MAX]);
void compute_temp (PROB_DATATYPE prob[NBETA_PER_WORD][NPROB_MAX], float *temp_beta_shared, int word);
void shuffle (int *temp_idx_shared, float *temp_beta_shared, int *E, int mod);
PROB_DATATYPE myrand ();

// kernel_l2.c
void mc (float *, Parameter, const int, const int);
void pt (int *, float *, Parameter, const int, const int);

// kernel_l1.c
void kernel_warmup (const Parameter);
void kernel_swap (const int, const Parameter);
void kernel_rearrange (const Parameter);
void kernel_compute_q (const int, const Parameter);



#endif /* SIM_H */

