#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "sim.h"
#include "include/yeah/timing.h"


void
host_launcher (float beta_low, float beta_high, char *mydir, int node, int device)
{
  // initilize random sequence
  srand (time (NULL) + 10 * node + device);


  // lattice
  MSC_DATATYPE *lattice;
  size_t lattice_sz = sizeof (MSC_DATATYPE) * SZ_CUBE * NWORD * REALIZATION;
  lattice = (MSC_DATATYPE *) malloc (lattice_sz);
  host_init_lattice (lattice);

  // lattice1
  // spins be rearranged to reflect the temperature order
  MSC_DATATYPE *lattice1;
  size_t lattice1_sz = sizeof (MSC_DATATYPE) * SZ_CUBE * REALIZATION;
  lattice1 = (MSC_DATATYPE *) malloc (lattice1_sz);

  // temp - index and beta
  Temp *temp;
  size_t temp_sz = sizeof (Temp) * NBETA_MAX * REALIZATION;
  temp = (Temp *) malloc (temp_sz);
  host_init_temp (temp, beta_low, beta_high);

  // st - status records
  St *st;
  size_t st_sz = sizeof (St) * ITER_SWAP / ITER_SWAP_KERN;
  st = (St *) malloc (st_sz);


  Parameter para;
  para.lattice = lattice;
  para.lattice1 = lattice1; 
  para.temp = temp;
  para.st = st;





  double t[4][2], t2[3], t3 = 0; // timing information

  char message[STR_LENG];

  putchar ('\n');
  host_report_speed_title ();


  // warm up runs
  t2[0] = HostTimeNow ();
  for (int i = 0; i < ITER_WARMUP; i += ITER_WARMUP_KERN) {
    t[0][0] = HostTimeNow ();

    kernel_warmup (para);

    t[0][1] = HostTimeNow ();
    sprintf (message, "n%03d d%d warmup %8d/%08d", node, device, i, ITER_WARMUP);
    host_report_speed (t[0][0], t[0][1], ITER_WARMUP_KERN, message);
  }

  t2[1] = HostTimeNow ();

  // swap runs
  for (int i = 0; i < ITER_SWAP; i += ITER_SWAP_KERN) {
    t[1][0] = HostTimeNow ();

    kernel_swap (i / ITER_SWAP_KERN, para);

    t[1][1] = HostTimeNow ();
    t3 += t[1][1] - t[1][0];

    kernel_rearrange (para);
    kernel_compute_q (i / ITER_SWAP_KERN, para);

    t[2][1] = HostTimeNow ();

    sprintf (message, "n%03d d%d PT     %8d/%08d", node, device, i, ITER_SWAP);
    host_report_speed (t[1][0], t[1][1], ITER_SWAP_KERN, message);
  }
  t2[2] = HostTimeNow ();



#ifndef NO_OUTPUT
  host_save_st (st, mydir, node, device);
#endif



  // report overall speed
  putchar ('\n');
  sprintf (message, "n%03d d%d overall warmup          ", node, device);
  host_report_speed (t2[0], t2[1], ITER_WARMUP, message);
  sprintf (message, "n%03d d%d overall PT (no measure) ", node, device);
  host_report_speed (0, t3, ITER_SWAP, message);
  sprintf (message, "n%03d d%d overall PT              ", node, device);
  host_report_speed (t2[1], t2[2], ITER_SWAP, message);
  sprintf (message, "n%03d d%d overall simulation      ", node, device);
  host_report_speed (t2[0], t2[2], ITER_WARMUP + ITER_SWAP, message);
  putchar ('\n');


  free (lattice);
  free (lattice1);
  free (temp);
  free (st);
}

