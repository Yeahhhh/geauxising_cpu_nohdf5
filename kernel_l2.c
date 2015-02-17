#include "sim.h"


void
mc (float *temp_beta_shared, Parameter para, const int realization, const int iter)
{
  MSC_DATATYPE *lattice = para.lattice;

  /// temperature scratchpad
  PROB_DATATYPE temp_prob_shared[NBETA_PER_WORD][NPROB_MAX];
  init_temp (temp_prob_shared);

  /// lattice scratchpad
  // sizeof(int32_t) * 16 * 16 * 16 = 16 KB
  MSC_DATATYPE l[SZ_CUBE];


  for (int word = 0; word < NWORD; ++word) {
    int lattice_offset = SZ_CUBE * NWORD * realization + SZ_CUBE * word;

    // initilize temperature scratchpad
    compute_temp (temp_prob_shared, temp_beta_shared, word);

    // import lattice scratchpad
    for (int i = 0; i < SZ_CUBE; ++i)
      l[i] = lattice[lattice_offset + i];



    for (int i = 0; i < iter; ++i) {

      // two phases update
      for (int run = 0; run < 2; run++) {
	for (int z = 0; z < L; ++z) {
	  const int za = (z + L - 1) % L;
	  const int zb = (z + 1) % L;
	  for (int y = 0; y < L; ++y) {
	    const int ya = (y + L - 1) % L;
	    const int yb = (y + 1) % L;
	    for (int x = ((y & 1) ^ (z & 1) ^ run); x < L; x += 2) {
	      const int xa = (x + L - 1) % L;
	      const int xb = (x + 1) % L;

	      MSC_DATATYPE c = l[CUBEIDX (z, y, x)];	// center
	      MSC_DATATYPE n0 = l[CUBEIDX (z, y, xa)];	// left
	      MSC_DATATYPE n1 = l[CUBEIDX (z, y, xb)];	// right
	      MSC_DATATYPE n2 = l[CUBEIDX (z, ya, x)];	// up
	      MSC_DATATYPE n3 = l[CUBEIDX (z, yb, x)];	// down
	      MSC_DATATYPE n4 = l[CUBEIDX (za, y, x)];	// front
	      MSC_DATATYPE n5 = l[CUBEIDX (zb, y, x)];	// back
	      
	      n0 = MASK_A * ((c >> SHIFT_J0) & 1) ^ n0 ^ c;
	      n1 = MASK_A * ((c >> SHIFT_J1) & 1) ^ n1 ^ c;
	      n2 = MASK_A * ((c >> SHIFT_J2) & 1) ^ n2 ^ c;
	      n3 = MASK_A * ((c >> SHIFT_J3) & 1) ^ n3 ^ c;
	      n4 = MASK_A * ((c >> SHIFT_J4) & 1) ^ n4 ^ c;
	      n5 = MASK_A * ((c >> SHIFT_J5) & 1) ^ n5 ^ c;
	      
	      for (int s = 0; s < NBETA_PER_SEG; ++s) {
		MSC_DATATYPE e =
		  ((n0 >> s) & MASK_S) +
		  ((n1 >> s) & MASK_S) +
		  ((n2 >> s) & MASK_S) +
		  ((n3 >> s) & MASK_S) +
		  ((n4 >> s) & MASK_S) + ((n5 >> s) & MASK_S);
		e = (e << 1) + ((c >> s) & MASK_S);
		MSC_DATATYPE flip = 0;
		
		for (int seg_offset = 0; seg_offset < SHIFT_MAX; seg_offset += NBIT_PER_SEG) {
		  PROB_DATATYPE val = temp_prob_shared[seg_offset + s][(e >> seg_offset) & MASK_E];
		  flip |= ((MSC_DATATYPE) (myrand () < val) << seg_offset);	// myrand < val ? 1 : 0;
		}
		c ^= (flip << s);
	      } // s

	    l[CUBEIDX (z, y, x)] = c;
	    }			// x
	  }                     // y
	}                       // z
      }				// run
    }				// i

    // export lattice scratchpad
    for (int i = 0; i < SZ_CUBE; ++i)
      lattice[lattice_offset + i] = l[i];

  }				// word

}






void
pt (int *temp_idx_shared, float *temp_beta_shared, Parameter para, const int realization, const int mod)
{
  MSC_DATATYPE *lattice = para.lattice;

  int E_shared[NBETA_PER_WORD];
  int E[NBETA];
  int Eh[NBETA];


  /// lattice scratchpad
  // sizeof (int32_t) * 16 * 16 * 16 = 16 KB
  MSC_DATATYPE l[SZ_CUBE];



  for (int word = 0; word < NWORD; ++word) {
    int lattice_offset = SZ_CUBE * NWORD * realization + SZ_CUBE * word;

    // import lattice scratchpad
    for (int i = 0; i < SZ_CUBE; ++i)
      l[i] = lattice[lattice_offset + i];

    // reset partial status
    for (int b = 0; b < NBETA_PER_WORD; b++)
      E_shared[b] = 0;


    for (int z = 0; z < L; ++z) {
      const int za = (z + L - 1) % L;
      const int zb = (z + 1) % L;
      for (int y = 0; y < L; ++y) {
	const int ya = (y + L - 1) % L;
	const int yb = (y + 1) % L;
	for (int x = 0; x < L; ++x) {
	  const int xa = (x + L - 1) % L;
	  const int xb = (x + 1) % L;

	  MSC_DATATYPE c = l[CUBEIDX (z, y, x)];	// center
	  MSC_DATATYPE n0 = l[CUBEIDX (z, y, xa)];	// left
	  MSC_DATATYPE n1 = l[CUBEIDX (z, y, xb)];	// right
	  MSC_DATATYPE n2 = l[CUBEIDX (z, ya, x)];	// up
	  MSC_DATATYPE n3 = l[CUBEIDX (z, yb, x)];	// down
	  MSC_DATATYPE n4 = l[CUBEIDX (za, y, x)];	// front
	  MSC_DATATYPE n5 = l[CUBEIDX (zb, y, x)];	// back
	  
	  n0 = MASK_A * ((c >> SHIFT_J0) & 1) ^ n0 ^ c;
	  n1 = MASK_A * ((c >> SHIFT_J1) & 1) ^ n1 ^ c;
	  n2 = MASK_A * ((c >> SHIFT_J2) & 1) ^ n2 ^ c;
	  n3 = MASK_A * ((c >> SHIFT_J3) & 1) ^ n3 ^ c;
	  n4 = MASK_A * ((c >> SHIFT_J4) & 1) ^ n4 ^ c;
	  n5 = MASK_A * ((c >> SHIFT_J5) & 1) ^ n5 ^ c;

	  for (int s = 0; s < NBETA_PER_SEG; s++) {
	    MSC_DATATYPE e =
	      ((n0 >> s) & MASK_S) +
	      ((n1 >> s) & MASK_S) +
	      ((n2 >> s) & MASK_S) +
	      ((n3 >> s) & MASK_S) +
	      ((n4 >> s) & MASK_S) +
	      ((n5 >> s) & MASK_S);

	    for (int seg_offset = 0; seg_offset < SHIFT_MAX; seg_offset += NBIT_PER_SEG) {
	      E_shared[seg_offset + s] += (e >> seg_offset) & MASK_E;	// range: [0,6]
	    }
	  }

	}			// x
      }				// y
    }				// z

    for (int b = 0; b < NBETA_PER_WORD; ++b)
      E[NBETA_PER_WORD * word + b] = E_shared[b];


    /// energy contribute by external field

    for (int b = 0; b < NBETA_PER_WORD; b++)
      E_shared[b] = 0;

    for (int z = 0; z < L; ++z) {
      for (int y = 0; y < L; ++y) {
	for (int x = 0; x < L; ++x) {
	  MSC_DATATYPE c = l[CUBEIDX (z, y, x)];

	  for (int i = 0; i < NSEG_PER_WORD; ++i) {
	    for (int j = 0; j < NBETA_PER_SEG; ++j) {
	      const int position = NBIT_PER_SEG * i + j;
	      E_shared[position] += ((c >> position) & 1);
	    } // j
	  } // i
	} // x
      } // y
    } // z

    
    for (int b = 0; b < NBETA_PER_WORD; ++b)
      Eh[NBETA_PER_WORD * word + b] = E_shared[b];

  }				// word;



  // convert E from [0,6] to [-6,6], e = e * 2 - 6
  // E = sum_CUBE (e * 2 - 6)
  //   = 2 * sum_CUBE e - 6 * SZ_CUBE

  // conver Eh from [0,1] to [-1,1], e = e * 2 - 1
  // Eh = 2 * sum_CUBE e - SZ_CUBE

  for (int b = 0; b < NBETA; ++b) {
    E[b] = E[b] * 2 - 6 * SZ_CUBE;
    Eh[b] = Eh[b] * 2 - SZ_CUBE;
    E[b] = E[b] + Eh[b] * H;
  }

  shuffle (temp_idx_shared, temp_beta_shared, E, mod);
}

