#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <math.h>

#include "sim.h"

#ifndef NO_OUTPUT
#include "hdf5.h"
#include <mpi.h>
#endif




/*
// convert lattice allocation, shared -> seperate
for (int run = 0; run < 2; run++) {
  for (z = 0; z < L; z++) {
    for (y = 0; y < L; y++) {
      for (x = 0; x < L_HF; x++) {
	idx_src = (L * L_HF * z) + (L_HF * y) + (x * 2 + run);
	idx_dst = (L * L_HF * z) + (L_HF * y) + (x + SZ_CUBE_HF * run);
	l1[idx_dst] = l0[idx_src];
      }
    }   
  }   
 }
*/




// initiate J
void
host_init_J (MSC_DATATYPE * l)
{
  MSC_DATATYPE jx[SZ_CUBE];
  MSC_DATATYPE jy[SZ_CUBE];
  MSC_DATATYPE jz[SZ_CUBE];

  MSC_DATATYPE j;
  int i;
  int x, xa, xb, y, ya, yb, z, za, zb;
  int idx;

  MSC_DATATYPE mask_s = MASK_S0;

  for (i = 0; i < SZ_CUBE; i++) {
    jx[i] = 1;
    jy[i] = 1;
    jz[i] = 1;
#ifdef RANDJ
    jx[i] = rand () & 1;
    jy[i] = rand () & 1;
    jz[i] = rand () & 1;
#endif
  }

  for (z = 0; z < L; z++) {
    za = (z + L - 1) % L;
    zb = z;			// zb != z + 1
    for (y = 0; y < L; y++) {
      ya = (y + L - 1) % L;
      yb = y;
      for (x = 0; x < L; x++) {
	xa = (x + L - 1) % L;
	xb = x;
	idx = L * L * z + L * y + x;

	j =			// integrated
	  MASK_J0 * jx[L * L * z + L * y + xa] |	// left
	  MASK_J1 * jx[L * L * z + L * y + xb] |	// right
	  MASK_J2 * jy[L * L * z + L * ya + x] |	// up
	  MASK_J3 * jy[L * L * z + L * yb + x] |	// down
	  MASK_J4 * jz[L * L * za + L * y + x] |	// front
	  MASK_J5 * jz[L * L * zb + L * y + x];	// back
	l[idx] = j | (l[idx] & mask_s);
      }
    }
  }

}


void
host_init_S (MSC_DATATYPE * l)
{

  MSC_DATATYPE s;
  MSC_DATATYPE mask_s = MASK_S0;


  for (int i = 0; i < SZ_CUBE; i++) {
    s = 1;
#ifdef RANDS
    s = rand ();
#endif
    s = s & mask_s;
    l[i] = (l[i] & MASK_J) | s;
  }
}



// logically, lattice[REALIZATION][SZ_CUBE]
void
host_init_lattice (MSC_DATATYPE * lattice)
{
  MSC_DATATYPE *l;		// local lattice
  l = (MSC_DATATYPE *) malloc (sizeof (MSC_DATATYPE) * SZ_CUBE);
  for (int i = 0; i < SZ_CUBE; i++) {
    l[i] = 0;
  }

  for (int cube = 0; cube < REALIZATION; cube++) {
    host_init_S (l);
    if ((cube & 1) == 0)	// common J for adjacent realization
      host_init_J (l);

    for (int word = 0; word < NWORD; word++) {
      int offset = (NWORD * SZ_CUBE * cube) + (SZ_CUBE * word);
      for (int i = 0; i < SZ_CUBE; i++)
	lattice[offset + i] = l[i];
    }
  }

  free (l);
}



// logically, lattice[REALIZATION][SZ_CUBE]
void
host_save_lattice (MSC_DATATYPE * lattice, char *mydir, int node, int device)
{
  for (int cube = 0; cube < REALIZATION; cube++) {
    char myfile[STR_LENG];
    sprintf (myfile, "%s/lattice_node%03d_dev%d_%04d.txt", mydir, node,
	     device, cube);
    FILE *fp = fopen (myfile, "wb");
    if (fp == NULL) {
      fprintf (stderr, "failed to open %s\n", myfile);
      exit (1);
    }
    fwrite (lattice, sizeof (MSC_DATATYPE), REALIZATION * SZ_CUBE * NWORD, fp);

    fclose (fp);
  }

}



// logically, lattice[REALIZATION][SZ_CUBE]
void
host_read_lattice (MSC_DATATYPE * lattice, char *mydir, int node, int device)
{
  for (int cube = 0; cube < REALIZATION; cube++) {
    char myfile[STR_LENG];
    sprintf (myfile, "%s/lattice_node%03d_dev%d_%04d.txt", mydir, node,
	     device, cube);
    FILE *fp = fopen (myfile, "rb");
    if (fp == NULL) {
      fprintf (stderr, "failed to open %s\n", myfile);
      exit (1);
    }
    fread (lattice, sizeof (MSC_DATATYPE), REALIZATION * SZ_CUBE * NWORD, fp);

    fclose (fp);
  }

}





/*
  output format of
  e[REALIZATION].txt
  q[REALIZATION_PAIR].txt

  X axis: beta, from low to high
  Y axis: iteration space
*/


void
host_save_st (St * st, char *mydir, int node, int device)      //, int nrank_per_node)
{

#ifndef NO_OUTPUT

  MPI_Group orig_group, new_group;
  MPI_Comm_group (MPI_COMM_WORLD, &orig_group);

  MPI_Comm new_comm;
  MPI_Info info = MPI_INFO_NULL;
  int mpi_rank;

  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
  MPI_Group_incl (orig_group, 1, &mpi_rank, &new_group);
  MPI_Comm_create (MPI_COMM_WORLD, new_group, &new_comm);

  hid_t file_id, dset_id, group_id;     /* file and dataset
                                           identifiers */
  hid_t filespace;              /* file and memory dataspace identifiers */
  hsize_t dimsf[] = { ITER_SWAP / ITER_SWAP_KERN, NBETA };  /* dataset dimensions */
  hid_t plist_id;               /* property list identifier */
  herr_t status;

  char myfile[STR_LENG];
  char name[500];


  float q[REALIZATION_HF][19][ITER_SWAP / ITER_SWAP_KERN][NBETA];



  for (int rec = 0; rec < ITER_SWAP / ITER_SWAP_KERN; rec++) {
    for (int pair = 0; pair < REALIZATION_HF; pair++) {
      for (int beta = 0; beta < NBETA; beta++) {
        q[pair][0][rec][beta] = st[rec].q[pair][beta];
        for (int m = 0; m < 3; m++) {
          q[pair][m + 1][rec][beta] = st[rec].qk_real[3][pair][beta];
          q[pair][m + 10][rec][beta] = st[rec].qk_imag[m][m][beta];
        }
	for (int m = 0; m < 6; m++) {
          q[pair][m + 4][rec][beta] = st[rec].qk2_real[m][pair][beta];
          q[pair][m + 13][rec][beta] = st[rec].qk2_imag[m][pair][beta];
        }
      }
    }
  }



  for (int pair = 0; pair < REALIZATION / 2; pair++) {

    sprintf (myfile, "%s/sample_node%03d_dev%d_%04d.h5", mydir, node, device,
             pair);
    plist_id = H5Pcreate (H5P_FILE_ACCESS);
    H5Pset_fapl_mpio (plist_id, new_comm, info);

    file_id = H5Fcreate (myfile, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose (plist_id);

    plist_id = H5Pcreate (H5P_DATASET_XFER);

    H5Pset_dxpl_mpio (plist_id, H5FD_MPIO_COLLECTIVE);
    //H5Pset_dxpl_mpio (plist_id, H5FD_MPIO_INDEPENDENT);

    filespace = H5Screate_simple (MPIRANK_PER_NODE, dimsf, NULL);
    sprintf (name, "q");
    dset_id = H5Dcreate (file_id, name, H5T_NATIVE_FLOAT, filespace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dclose (dset_id);
    H5Sclose (filespace);

    dset_id = H5Dopen (file_id, name, H5P_DEFAULT);
    status =
      H5Dwrite (dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, plist_id,
                q[pair][0]);
    H5Dclose (dset_id);

    for (int m = 1; m < 10; m++) {
      filespace = H5Screate_simple (MPIRANK_PER_NODE, dimsf, NULL);
      sprintf (name, "qk_real_%02d", m);
      dset_id = H5Dcreate (file_id, name, H5T_NATIVE_FLOAT, filespace,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dclose (dset_id);
      H5Sclose (filespace);

      dset_id = H5Dopen (file_id, name, H5P_DEFAULT);
      status =
        H5Dwrite (dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, plist_id,
                  q[pair][m]);
      H5Dclose (dset_id);

    }


    for (int m = 1; m < 10; m++) {
      filespace = H5Screate_simple (MPIRANK_PER_NODE, dimsf, NULL);

      sprintf (name, "qk_imag_%02d", m);
      dset_id = H5Dcreate (file_id, name, H5T_NATIVE_FLOAT, filespace,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dclose (dset_id);
      H5Sclose (filespace);

      dset_id = H5Dopen (file_id, name, H5P_DEFAULT);
      status =
        H5Dwrite (dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, plist_id,
                  q[pair][m + 9]);
      H5Dclose (dset_id);
    }


    H5Pclose (plist_id);
    H5Fclose (file_id);

  }
#endif
}





#if 0
void
host_save_st (St * st, char *mydir, int node, int device)
{
  //printf ("node%03d device%d write outputs to %s/\n", node, device, mydir);

  const int rec_max = ITER_SWAP / ITER_SWAP_KERN;

  /*
     for (int realization = 0; realization < REALIZATION; realization++) {
     char myfile[STR_LENG];
     sprintf (myfile, "%s/node%03d_dev%d_e%04d.txt", mydir, node, device, realization);
     FILE *fp = fopen (myfile, "a");
     if (fp == NULL) {
     fprintf (stderr, "failed to open %s\n", myfile);
     exit (1);
     }

     for (int rec = 0; rec < rec_max; rec++) {
     for (int b = 0; b < NBETA; b++)
     fprintf (fp, "%8.2f ", st[rec].e[realization][b]);
     fprintf (fp, "\n");
     }

     fclose (fp);

     }
   */

  //output q(0)
  for (int pair = 0; pair < REALIZATION / 2; pair++) {
    char myfile[STR_LENG];
    sprintf (myfile, "%s/q_node%03d_dev%d_%04d.txt", mydir, node, device,
	     pair);
    FILE *fp = fopen (myfile, "a");
    if (fp == NULL) {
      fprintf (stderr, "failed to open %s\n", myfile);
      exit (1);
    }

    for (int rec = 0; rec < rec_max; rec++) {
      for (int b = 0; b < NBETA; b++)
	fprintf (fp, "%10.3f", st[rec].q[pair][b]);
      fprintf (fp, "\n");
    }

    fclose (fp);
  }

  for (int j = 0; j < 3; j++) {
    //output q(k1)
    for (int pair = 0; pair < REALIZATION / 2; pair++) {

      char myfile[STR_LENG];
      sprintf (myfile, "%s/qk1r_%d_node%03d_dev%d_%04d.txt", mydir, j, node,
	       device, pair);
      FILE *fp = fopen (myfile, "a");
      if (fp == NULL) {
	fprintf (stderr, "failed to open %s\n", myfile);
	exit (1);
      }

      for (int rec = 0; rec < rec_max; rec++) {
	for (int b = 0; b < NBETA; b++)
	  fprintf (fp, "%10.3f", st[rec].qk_real[j][pair][b]);
	fprintf (fp, "\n");
      }

      fclose (fp);
    }

    for (int pair = 0; pair < REALIZATION / 2; pair++) {
      char myfile[STR_LENG];
      sprintf (myfile, "%s/qk1i_%d_node%03d_dev%d_%04d.txt", mydir, j, node,
	       device, pair);
      FILE *fp = fopen (myfile, "a");
      if (fp == NULL) {
	fprintf (stderr, "failed to open %s\n", myfile);
	exit (1);
      }

      for (int rec = 0; rec < rec_max; rec++) {
	for (int b = 0; b < NBETA; b++)
	  fprintf (fp, "%10.3f", st[rec].qk_imag[j][pair][b]);
	fprintf (fp, "\n");
      }

      fclose (fp);
    }
  }

  for (int j = 0; j < 6; j++) {
    //output q(k2)
    for (int pair = 0; pair < REALIZATION / 2; pair++) {
      char myfile[STR_LENG];
      sprintf (myfile, "%s/qk2r_%d_node%03d_dev%d_%04d.txt", mydir,j ,node, device,
	       pair);
      FILE *fp = fopen (myfile, "a");
      if (fp == NULL) {
	fprintf (stderr, "failed to open %s\n", myfile);
	exit (1);
      }

      for (int rec = 0; rec < rec_max; rec++) {
	for (int b = 0; b < NBETA; b++)
	  fprintf (fp, "%10.3f", st[rec].qk2_real[j][pair][b]);
	fprintf (fp, "\n");
      }

      fclose (fp);
    }
    for (int pair = 0; pair < REALIZATION / 2; pair++) {
      char myfile[STR_LENG];
      sprintf (myfile, "%s/qk2i_%d_node%03d_dev%d_%04d.txt", mydir, j,node, device,
	       pair);
      FILE *fp = fopen (myfile, "a");
      if (fp == NULL) {
	fprintf (stderr, "failed to open %s\n", myfile);
	exit (1);
      }

      for (int rec = 0; rec < rec_max; rec++) {
	for (int b = 0; b < NBETA; b++)
	  fprintf (fp, "%10.3f", st[rec].qk2_imag[j][pair][b]);
	fprintf (fp, "\n");
      }

      fclose (fp);
    }
  }
}
#endif






// initialize temps
void
host_init_temp (Temp * temp, float beta_low, float beta_high)
{
  //const float beta_delta = (beta_high - beta_low) / NBETA;
  double a = beta_low;
  Temp *tmp = (Temp *) malloc (sizeof (Temp) * NBETA_MAX);
  const double beta_ratio = exp (log (beta_high / beta_low) / (NBETA - 1));

  //printf ("List of Betas used in this run:\n");
  // generating temperatures for one lattice
  for (int b = 0; b < NBETA; b++) {
    tmp[b].idx = b;
    tmp[b].beta = a;
    //a += beta_delta;
    a *= beta_ratio;
    //printf ("%f\n", tmp[b].beta);
  }

  for (int b = NBETA; b < NBETA_MAX; b++) {
    tmp[b].idx = b;
    tmp[b].beta = 0;
  }

  // duplication for all lattices
  for (int lattice = 0; lattice < REALIZATION; lattice++) {
    for (int b = 0; b < NBETA_MAX; b++) {
      int i = NBETA_MAX * lattice + b;
      //temp_idx[i] = tmp[b].idx;
      //temp_beta[i] = tmp[b].beta;
      temp[i] = tmp[b];
    }
  }

  free (tmp);
}




void
host_report_speed_title ()
{
  printf ("\t\t\t\t     second   iter/s    spin/s    s/spin    ");
  printf ("ps/spin   ");
  //printf ("SM_R/W    ");
  //printf ("GFLOPS  ");
  putchar ('\n');
}



void
host_report_speed (double start, double stop, int iter, char *event)
{
  double s = stop - start;
  float iter_per_s = iter / s;
  float spin_per_s = NBETA * SZ_CUBE * REALIZATION * iter_per_s;
  float s_per_spin = 1 / spin_per_s;
  //int sm_w_bw = NWORD * sizeof (MSC_DATATYPE) * SZ_CUBE * iter_per_s / 1048576;
  //int sm_r_bw = sm_w_bw * 7;  // average per TB, GB/s

  printf ("%s   %9.6f, ", event, s);
  printf ("%.3e %.3e %.3e ", iter_per_s, spin_per_s, s_per_spin);
  printf ("%06.3f    ", s_per_spin * 1.0e12);
  //printf ("%04d/%04d ", sm_r_bw, sm_w_bw);
  //printf ("%s", "???");
  putchar ('\n');
}





void
host_usage (char *bin)
{
  fprintf (stderr, "usage: %s [options]\n", bin);
  fprintf (stderr, " -l <beta_lower_bound>     default %f\n", BETA_LOW);
  fprintf (stderr, " -u <beta_upper_bound>     default %f\n", BETA_HIGH);
  fprintf (stderr, " where lower bound must less than upper bound\n");
  exit (1);
}



void
host_summary (float beta_low, float beta_high, char *mydir)
{
  printf ("\n");
  printf ("lattice size:\t\t %d*%d*%d\n", L, L, L);
  printf ("realizations per GPU:\t %d\n", REALIZATION);
  printf ("parallel betas:\t\t %.5f-%.5f,", beta_low, beta_high);
  printf (" %d samples\n", NBETA);
  printf ("external field:\t\t %.5f\n", H);
  printf ("ITER_WARMUP:\t\t %d\n", ITER_WARMUP);
  printf ("ITER_WARMUP_KERN:\t %d\n", ITER_WARMUP_KERN);
  printf ("ITER_WARMUP_KERNFUNC:\t %d\n", ITER_WARMUP_KERNFUNC);
  printf ("ITER_SWAP:\t\t %d\n", ITER_SWAP);
  printf ("ITER_SWAP_KERN:\t\t %d\n", ITER_SWAP_KERN);
  printf ("ITER_SWAP_KERNFUNC:\t %d\n", ITER_SWAP_KERNFUNC);
  printf ("output dir:\t\t %s\n", mydir);
  printf ("\n");

//  printf ("threads per block:\t %d\n", BD);
//  printf ("blocks per GPU:\t\t %d\n", REALIZATION);
  printf ("\n");
}

