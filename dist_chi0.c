#include <stdio.h>
#include <stdlib.h>


#include "hdf5.h"
#include "string.h"

#include "sim.h"
/*
#define SZ 16
#define SZ_CUBE (SZ*SZ*SZ)
#define ITER_SWAP 20000000
#define ITER_SWAP_KERN 1000
#define NBETA 24
#define STR_LENG 64
*/


// analysis realization_[REALIZATION]_e.txt

//#define RUNS 60
//#define SAMPLES (RUNS*REALIZATION/2)
//#define SAMPLES 6912

#define rec_max ( ITER_SWAP / ITER_SWAP_KERN)
//#define SAMPLES 1600

#define RANKC        1



int
main (int argc, char **argv)
{
  char myh5[STR_LENG] = { 0 };

  int rec_start = rec_max/2;
  float *qr;
  float k1r,x1,k1r2,k1r3,k1r4;
  char name[500];

  hid_t memspace;
  hsize_t count[2];
  hsize_t offset[2];
  hsize_t col_dims[1];
  hid_t file_id, dataset_id, dataspace_id;
  hid_t filespace;
  herr_t status;

  offset[0] = 0;

  count[0] = rec_max;
  count[1] = 1;
  col_dims[0] = rec_max;

  qr = (float *) malloc (sizeof (float) * rec_max);


  int pair, rec, k, g;

  if (argc == 3)
    {
      pair=atoi(argv[1]); 
      offset[1]=atoi(argv[2]);
    }
  else
    {
      printf ("Error: usage: dist_chi0 sample_index beta_index\n");
      return 0;
    }



  //  for (pair = 0; pair < SAMPLES; pair++) {

    sprintf(myh5,"sample_pair_%05d.h5",pair);

    file_id = H5Fopen (myh5, H5F_ACC_RDONLY, H5P_DEFAULT);

    if(file_id<0){
      printf("File do not exist!\n");
      exit(0);
    }



    x1 = 0.0f;
    k1r = 0.0f;
    k1r2 = 0.0f;
    k1r3 = 0.0f;
    k1r4 = 0.0f;

    sprintf (name, "q");      
    dataset_id = H5Dopen2 (file_id, name, H5P_DEFAULT);
    
    filespace = H5Dget_space (dataset_id);
    memspace = H5Screate_simple (RANKC, col_dims, NULL);
    status = H5Sselect_hyperslab (filespace, H5S_SELECT_SET, offset, NULL,
				    count, NULL);
    status =
      H5Dread (dataset_id, H5T_NATIVE_FLOAT, memspace, filespace,
	       H5P_DEFAULT, qr);
    status = H5Dclose (dataset_id);


  

    for (rec = rec_start; rec < rec_max; rec++) {
      double q=qr[rec]/SZ_CUBE;
      k1r += q;
      k1r2 += q*q;
      k1r3 += q*q*q;
      k1r4 += q*q*q*q;
      
    }
    
    x1=k1r2;

    double ax1 = x1 / (rec_max - rec_start);
    double aq1r = k1r / (rec_max - rec_start);
    double aq1r2 = k1r2 / (rec_max - rec_start);
    double aq1r3 = k1r3 / (rec_max - rec_start);
    double aq1r4 = k1r4 / (rec_max - rec_start);

    ax1 = (ax1 - aq1r * aq1r)*SZ_CUBE;

    status = H5Fclose (file_id);

    printf ("%7d\t%2d\t%1.5f\t%1.5f\t%1.5f\t%1.5f\t%1.5f\n", pair,(int)offset[1], ax1,aq1r,aq1r2,aq1r3,aq1r4);


  return 0;
}
