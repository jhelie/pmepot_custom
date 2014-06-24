/*
 * Top level PME potential routines
 *
 * $Id: pmepot.c,v 1.4 2005/07/20 15:37:39 johns Exp $
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

struct pmepot_data_struct {
  int dims[5];
  int grid_size;
  int max_dim;
  int fft_ntable;
  float ewald_factor;
  float oddd[12];
  int avg_count;
  float *avg_potential;
  float *fft_table;
  float *fft_work;
};

#include "pmepot.h"
#include "pub3dfft.h"

pmepot_data* pmepot_create(int *dims, float ewald_factor) {
  //called by:  set pme [pmepot_create $grid $ewaldfactor]
  
  pmepot_data *data;
  int grid_size, max_dim;

  if ( dims[0] < 8 ) return 0;
  if ( dims[1] < 8 ) return 0;
  if ( dims[2] < 8 ) return 0;
  if ( dims[2] % 2 ) return 0;
  if ( ewald_factor <= 0. ) return 0;

  data = malloc(sizeof(pmepot_data));
  if ( ! data ) return 0;

  //dims= [gridA gridB gridC]
  data->avg_count = 0;			//no frame processed yet
  data->ewald_factor = ewald_factor;
  data->dims[0] = dims[0];  	//gridA
  data->dims[1] = dims[1];		//gridB
  data->dims[2] = dims[2];		//gridC
  data->dims[3] = dims[1];		//gridB
  data->dims[4] = dims[2] + 2;	//gridC+2 : gives himself a bit of room in the z direction?
  grid_size = data->dims[0] * data->dims[3] * data->dims[4];	//gridA*gridB*(gridC+2)
  data->grid_size = grid_size;	//number of grid points
  max_dim = dims[0] > dims[1] ? dims[0] : dims[1];
  max_dim = max_dim > dims[2] ? max_dim : dims[2];
  data->max_dim = max_dim;
  data->fft_ntable = 4*max_dim+15;

  //allocates enough memory to store a grid of floats: avg_potential is the float values on the grid
  data->avg_potential = malloc(grid_size * sizeof(float));
  
  data->fft_table = malloc(3 * data->fft_ntable * sizeof(float));
  data->fft_work = malloc(2 * max_dim * sizeof(float));
  if ( ! data->avg_potential || ! data->fft_table || ! data->fft_work ) {
    if ( data->avg_potential) free(data->avg_potential);
    if ( data->fft_table) free(data->fft_table);
    if ( data->fft_work) free(data->fft_work);
    free(data);
    return 0;
  }

  //call a function from pub3dfft.c
  pubd3di(dims[2], dims[1], dims[0], data->fft_table, data->fft_ntable);

  //return the object data created
  return data;
}

void pmepot_destroy(pmepot_data *data) {
  free(data->avg_potential);
  free(data->fft_table);
  free(data->fft_work);
  free(data);
}

void scale_grid(const int *dims, float *arr, const float factor) {
  int grid_size, i;
  grid_size = dims[0] * dims[3] * dims[4];
  for ( i=0; i<grid_size; ++i ) {
	  arr[i] *= factor;
	  //printf("-index: %d\n", i);
	  //printf("-potential: %f\n", arr[i]);
	  }
}

void add_to_grid(const int *dims, const int count,float *avg, const float *arr) {
  
  //called by: add_to_grid(data->dims,data->avg_potential,q_arr);
  //so this is doing:
  //	data->avg_potential[i]+=q_arr[i] for all grid points
  
  int grid_size, i;
  grid_size = dims[0] * dims[3] * dims[4];
  
  if (count==0) {
	for ( i=0; i<grid_size; ++i ) {
		avg[i]=0;
	    //printf("-index: %d\n", i);
		//printf("-potential before addition: %f\n", avg[i]);
		//printf("-potential to be added: %f\n", arr[i]);
		avg[i] += arr[i];
		//printf("-potential after addition: %f\n", avg[i]);
	}
  } else {
		for ( i=0; i<grid_size; ++i ) {
			//printf("-count: %d\n", count);
			//printf("-index: %d\n", i);
			//printf("-potential before addition: %f\n", avg[i]);
			//printf("-potential to be added: %f\n", arr[i]);
			avg[i] += arr[i];
			//printf("-potential after addition: %f\n", avg[i]);
		}
  }
}

void display_grid(const int *dims, float *grid) {
  
  //called by: add_to_grid(data->dims,data->avg_potential,q_arr);
  //so this is doing:
  //	data->avg_potential[i]+=q_arr[i] for all grid points

  float min, max;
  min=0;
  max=0;  
  printf("-1ere passe de display\n");
  int grid_size, i;
  grid_size = dims[0] * dims[3] * dims[4];
  for ( i=0; i<grid_size; ++i ) {
	if (grid[i]>max) {
		printf("-index: %d\n", i);
		printf("-max_tmp: %1.16f\n", grid[i]);
		max=grid[i];
	}
	if (grid[i]<min) {
		printf("-index: %d\n", i);
		printf("-min_tmp: %1.16f\n", grid[i]);
		min=grid[i];
	}	
	if (grid[i]>200) {
		printf("-index: %d\n", i);
		printf("- potential huge max: %1.16f\n", grid[i]);
		grid[i]=0;
	}
	if (grid[i]<-200) {
		printf("-index: %d\n", i);
		printf("-potential huge min: %1.16f\n", grid[i]);
		grid[i]=0;
	}

	if ( isnan(grid[i]) ) {
		printf("-index: %d\n", i);
		printf("- potential nan: %1.16f\n", grid[i]);
		grid[i]=0;
	}
  }
  printf("-max: %1.16f\n", max);
  printf("-min: %1.16f\n", min);

}

void display_grid2(const int *dims, float *grid) {
  
  //called by: add_to_grid(data->dims,data->avg_potential,q_arr);
  //so this is doing:
  //	data->avg_potential[i]+=q_arr[i] for all grid points
  
  float min2, max2;
  min2=0;
  max2=0;
  printf("-2eme passe de display\n");
  int grid_size, i;
  grid_size = dims[0] * dims[3] * dims[4];
  for ( i=0; i<grid_size; ++i ) {
	if (grid[i]>max2) {
		printf("-index: %d\n", i);
		printf("-max_tmp: %1.16f\n", grid[i]);
		max2=grid[i];
	}
	if (grid[i]<min2) {
		printf("-index: %d\n", i);
		printf("-min_tmp: %1.16f\n", grid[i]);
		min2=grid[i];
	}	
	if (grid[i] > 200 ) {
		printf("-index: %d\n", i);
		printf("-display huge: %1.16f\n", grid[i]);
	}
	if (grid[i] < -200 ) {
		printf("-index: %d\n", i);
		printf("-display huge: %1.16f\n", grid[i]);
	}

	if ( isnan(grid[i]) ) {
		printf("-index: %d\n", i);
		printf("-display nan: %1.16f\n", grid[i]);
	}
  }
  printf("-max2: %1.16f\n", max2);
  printf("-min2: %1.16f\n", min2);
}


int fill_charges(const int *dims, const float *cell, int natoms,	const float *xyzq, float *q_arr_charge, float *q_arr_recip, float *q_arr_direct, float *rcell, float *oddd);
//int fill_charges(const int *dims, const float *cell, int natoms,	const float *xyzq, float *q_arr_charge, float *q_arr_recip, float *rcell, float *oddd);

float compute_energy_recip(float *q_arr_recip, const float *cell, const float *rcell, const int *dims, float ewald);

int compute_energy_direct(float *q_arr_charge, float *q_arr_direct, const float *cell, const int *dims, float ewald);

#define COLOUMB 332.0636
#define BOLTZMAN 0.001987191

// This is the function which is run for each frame to calculate the potential
//****************************************************************************
int pmepot_add(pmepot_data *data, const float *cell, int natoms, const float *atoms) {
  
  //called in pmepot.tcl by:
  //pmepot_add $pme $cell [$sel get {x y z charge}]
  // data = pme
  // cell = cell
  // natoms and atoms =? [$sel get {x y z charge}]
  //explanation: actually called by tcl_pmepot_add in tcl_pmepot.c which performs some check and processess [$sel get {x y z charge}] to extract the natoms and atoms variables
  // natoms=number of atoms (int)
  // atoms= [natomsx4] float array (coords+charge for each atom)
  
  float *q_arr_recip;
  float *q_arr_direct;
  float *q_arr_charge;
  float rcell[12];

  //q_arr is a float array of charges (same size as avg_potential) NB: not sure what's the difference between the 2 if scenario?
  //q_arr is modified by the routines below so as to actually store not the charge but the potential values
  q_arr_charge = malloc(data->grid_size * sizeof(float));
  q_arr_recip = malloc(data->grid_size * sizeof(float));
  q_arr_direct = malloc(data->grid_size * sizeof(float));
  if ( ! q_arr_recip ) return -1;

  //calculate the interpolated charge for each point of the grid
  //************************************************************
  fill_charges(data->dims,cell,natoms,atoms,q_arr_charge,q_arr_recip,q_arr_direct,rcell,data->oddd);

  //calculate potential in reciprocal space
  //***************************************
  
  //does 3D FFT of the interpolated charge distribution on the grid: stores the coefficient in the same grid
  pubdz3d(1, data->dims[2], data->dims[1], data->dims[0], q_arr_recip, data->dims[4], data->dims[3], data->fft_table, data->fft_ntable, data->fft_work);

  //calculate coefficient of fourier series (i.e. potential in reciprocral space):
  //compute_energy_recip(q_arr_recip, cell, rcell, data->dims, data->ewald_factor);

  //does FFT inverse, i.e. calculate potential in direct space using the fourier coefficient
  //pubzd3d(-1, data->dims[2], data->dims[1], data->dims[0], q_arr_recip, data->dims[4], data->dims[3], data->fft_table, data->fft_ntable, data->fft_work);

  //convert into kT unit: the COLOUMB constant converts an energy in q^2.A^-1 into KCal.mol^-1 and 300*Boltzmann corresponds to the value of kT in KCal.mol^-1 at 300K 
  //scale_grid(data->dims,q_arr_recip,COLOUMB/(300.0*BOLTZMAN));
  //printf("scaling factor:%1.16f\n", COLOUMB/(300.0*BOLTZMAN));

  //calculate potential in direct space
  //***********************************
  
  //calculate direct potential for each grid_point
  compute_energy_direct(q_arr_charge, q_arr_direct, cell, data->dims, data->ewald_factor);
  
  //scale it
  //scale_grid(data->dims,q_arr_direct,COLOUMB/(300.0*BOLTZMAN));
  
  //add potentials to the grid
  //**************************
  //add_to_grid(data->dims,data->avg_count,data->avg_potential,q_arr_recip);
  add_to_grid(data->dims,data->avg_count,data->avg_potential,q_arr_direct);
  //add_to_grid(data->dims,data->avg_count,data->avg_potential,q_arr_charge);
  //display_grid(data->dims,data->avg_potential);
  //display_grid2(data->dims,data->avg_potential);
  free(q_arr_recip);
  free(q_arr_charge); 
  free(q_arr_direct); 
   
  //update the number of frames taken into account (needed for later averaging)
  data->avg_count += 1;

  return 0;
}

int write_dx_grid(FILE *file, const int *dims, const float *oddd, const float *data, float scale, const char *label);

int pmepot_writedx(pmepot_data *data, const char *filename) {
  FILE *file;
  int rval;
  if ( ! data->avg_count ) return -1;
  file = fopen(filename,"w");
  if ( ! file ) return -2;
  rval = write_dx_grid(file,data->dims,data->oddd,data->avg_potential,1./data->avg_count,"PME potential (kT/e, T=300K)");
  fclose(file);
  return rval * 10;
}

