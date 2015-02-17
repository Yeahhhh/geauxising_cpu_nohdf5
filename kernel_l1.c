#include "COPYING"

#include <math.h>
#include "sim.h"


void
kernel_warmup (const Parameter para)
{

  /// temperature
  // (4 * 32) * 2 = 256 B
  float temp_beta_shared[NBETA];
  for (int r = 0; r < REALIZATION; ++r) {
    for (int b = 0; b < NBETA; ++b) {
      temp_beta_shared[b] = para.temp[NBETA_MAX * r + b].beta;
    }
  }
  
  for (int r = 0; r < REALIZATION; ++r) {
    for (int i = 0; i < ITER_WARMUP_KERN; i += ITER_WARMUP_KERNFUNC) {
      mc (temp_beta_shared, para, r, ITER_WARMUP_KERNFUNC);
    }
  }
}




void
kernel_swap (const int rec, const Parameter para)
{
  Temp *temp = para.temp;

  /// temperature
  // (4 * 32) * 2 = 256 B
  int temp_idx_shared[NBETA];
  float temp_beta_shared[NBETA];

  // load temperature
  for (int r = 0; r < REALIZATION; ++r) {
    for (int b = 0; b < NBETA; ++b) {
      temp_idx_shared[b] = temp[NBETA_MAX * r + b].idx;
      temp_beta_shared[b] = temp[NBETA_MAX * r + b].beta;
    }
  }

  for (int r = 0; r < REALIZATION; ++r) {
    for (int i = 0; i < ITER_SWAP_KERN; i += ITER_SWAP_KERNFUNC) {
      int swap_mod = (i / ITER_SWAP_KERNFUNC) & 1;
      pt (temp_idx_shared, temp_beta_shared, para, r, swap_mod);
      mc (temp_beta_shared, para, r, ITER_SWAP_KERNFUNC);
    }
  }

  // store temperature
  for (int r = 0; r < REALIZATION; ++r) {
    for (int b = 0; b < NBETA; ++b) {
      temp[NBETA_MAX * r + b].idx = temp_idx_shared[b];
      temp[NBETA_MAX * r + b].beta = temp_beta_shared[b];
    }
  }

  // store energy status
  // for (int r = 0; r < REALIZATION; ++r)
  // for (int b = 0; b < NBETA; ++b)
  // para.st[rec].e[r][temp_idx_shared[b]] = E[b];
}






// rearrange the spins so that they matches the temperature order
// least significant bit - lowest temperature
// higher order bits - higher temperature
// for ACMSC

void
kernel_rearrange (const Parameter para)
{
  MSC_DATATYPE *lattice = para.lattice;
  MSC_DATATYPE *lattice1 = para.lattice1;
  Temp *temp = para.temp;

  // temperature scratchpad
  int temp_idx_shared[NBETA_PER_WORD];


  for (int r = 0; r < REALIZATION; ++r) {

    // initilize lattice1
    for (int site = 0; site < SZ_CUBE; ++site)
      lattice1[SZ_CUBE * r + site] = 0;

    const int word = 0;
    const int lattice_offset = (SZ_CUBE * NWORD * r) + (SZ_CUBE * word);

    for (int b = 0; b < NBETA_PER_WORD; ++b) {
      temp_idx_shared[b] = temp[NBETA_MAX * r + NBETA_PER_WORD * word + b].idx;

      for (int site = 0; site < SZ_CUBE; ++site) {
	MSC_DATATYPE oldword = lattice[lattice_offset + site];
	MSC_DATATYPE newword = 0;
	
	for (int i = 0; i < NSEG_PER_WORD; ++i) {
	  for (int j = 0; j < NBETA_PER_SEG; ++j) {
	    const int position = NBIT_PER_SEG * i + j;
	    const int b = NBETA_PER_SEG * i + j;
	    MSC_DATATYPE tmp = oldword >> position & 1;
	    tmp <<= temp_idx_shared[b];
	    newword |= tmp;
	  }
	}
	lattice1[SZ_CUBE * r + site] |= newword;
      }
    }
  }

}








void
kernel_compute_q (const int rec, const Parameter para)
{
  MSC_DATATYPE *lattice1 = para.lattice1;

  MSC_DATATYPE l1[SZ_CUBE];
  double qk_real[3][NBETA];
  double qk_imag[3][NBETA];
  double qk2_real[6][NBETA];
  double qk2_imag[6][NBETA];


  for (int r = 0; r < REALIZATION; ++r) {

    const int lattice_offset0 = SZ_CUBE * (r << 1);
    const int lattice_offset1 = lattice_offset0 + SZ_CUBE;
    const double k = 2 * PI / L;

    for (int site = 0; site < SZ_CUBE; ++site) {
      l1[site] =		// xord_word 
	lattice1[lattice_offset0 + site] ^
	lattice1[lattice_offset1 + site];
    }

    for (int b = 0; b < NBETA; ++b) {
      float q0 = 0.0f;
      for (int j = 0; j < 3; j++) {
	qk_real[j][b] = 0.0f;
	qk_imag[j][b] = 0.0f;
      }
      for (int j = 0; j < 6; j++) {
	qk2_real[j][b] = 0.0f;
	qk2_imag[j][b] = 0.0f;
      }

      MSC_DATATYPE xor_word;
      int xor_bit;

      for (int i = 0; i < SZ_CUBE; i++) {
	xor_word = l1[i];
	xor_bit = (xor_word >> b) & 0x1;
	xor_bit = 1 - (xor_bit << 1);	// parallel: +1, reverse: -1

	double bit = xor_bit;
	double x = i % L;
	double y = (i / L) % L;
	double z = (i / L) / L;
	/*      // 2 * pi / L * x_i
		angel1 = (double) (i % L) * 2 * PI / L;
		// 2 * pi / L * (x_i + y_i)
		angel2 = (double) (i % L + (i / L) % L) * 2 * PI / L;
	*/
	q0 += bit;
	/*
	  qk_real += (float)xor_bit * cos (angel1);
	  qk_imag += (float)xor_bit * sin (angel1);
	  qk2_real += (float)xor_bit * cos (angel2);
	  qk2_imag += (float)xor_bit * sin (angel2);
	*/
	qk_real[0][b] += bit * cos (x * k);
	qk_real[1][b] += bit * cos (y * k);
	qk_real[2][b] += bit * cos (z * k);
	
	qk_imag[0][b] += bit * sin (x * k);
	qk_imag[1][b] += bit * sin (y * k);
	qk_imag[2][b] += bit * sin (z * k);
	
	qk2_real[0][b] += bit * cos (x * k + y * k);
	qk2_real[1][b] += bit * cos (x * k - y * k);
	qk2_real[2][b] += bit * cos (x * k + z * k);
	qk2_real[3][b] += bit * cos (x * k - z * k);
	qk2_real[4][b] += bit * cos (y * k + z * k);
	qk2_real[5][b] += bit * cos (y * k - z * k);
	
	qk2_imag[0][b] += bit * sin (x * k + y * k);
	qk2_imag[1][b] += bit * sin (x * k - y * k);
	qk2_imag[2][b] += bit * sin (x * k + z * k);
	qk2_imag[3][b] += bit * sin (x * k - z * k);
	qk2_imag[4][b] += bit * sin (y * k + z * k);
	qk2_imag[5][b] += bit * sin (y * k - z * k);
      }


      St *st = para.st;

      // save measurements in "st"
      st[rec].q[r][b] = q0;
      for (int j = 0; j < 3; j++) {
	st[rec].qk_real[j][r][b] = qk_real[j][b];
	st[rec].qk_imag[j][r][b] = qk_imag[j][b];
      }
      //      (float) sqrt (qk_real * qk_real + qk_imag * qk_imag);
      for (int j = 0; j < 6; j++) {
	st[rec].qk2_real[j][r][b] = qk2_real[j][b];
	st[rec].qk2_imag[j][r][b] = qk2_imag[j][b];
      }
      //      (Float) sqrt (qk2_real * qk2_real + qk2_imag * qk2_imag);


    } // beta
  } // r
}
