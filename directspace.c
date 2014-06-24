//author: Jean Helie
//date: March 2013

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#if !defined(M_PI)
#define M_PI 3.14159265358979323846
#endif

//This function is called by the pmepot_add function for each frame
//*****************************************************************
float compute_energy_direct(float *q_arr_charge, float *q_arr_direct, const float *cell, const int *dims, float ewald) {
 
  //declare variables 
  int K1, K2, K3, dim2, dim3;
  int k1, k1_s, k1_sn, k2, k2_s, k2_sn, k3, k3_s, k3_sn;
  int i, j, l;
  int ind_current, ind_neighbour, check_origin;
  float C, i_4pi_e0;
  float dist_neighbour;
  int counter;
  
  //stores number of grid points in each direction in K1 (gridA), K2(gridB) and K3(gridC)
  K1=dims[0]; K2=dims[1]; K3=dims[2]; dim2=dims[3]; dim3=dims[4];

  //stores the value of the term 1/(4*Pi*e0)
  i_4pi_e0=0.898775518; //V.A.C-1;

  //case where the cell is a cuboid. i.e. A, B and C vectors are all orthogonal to each other (other cases not dealt with for now)
  if ( cell[4] == 0. && cell[5] == 0. && cell[6] == 0. && cell[8] == 0. && cell[9] == 0. && cell[10] == 0. ) {

	//debug
	/*printf("-grid x: %1.16f",cell[3]/(K1-1));
	printf("-grid y: %1.16f",cell[7]/(K2-1));
	printf("-grid z: %1.16f",cell[11]/(K3-1));*/
	
	//browse all the grid points
    counter=0;
    for ( k1=0; k1<K1; ++k1 ) {
		for ( k2=0; k2<K2; ++k2 ) {
			for ( k3=0; k3<K3; ++k3 ) {
				//calculate index number of current grid point
				ind_current=k1*dim2*dim3+k2*dim3+k3;
				//printf("-ind_current: %d\n", ind_current);
				
				//shift grid coordinates to ease loop browsing below
				k1_s = k1-8;
				k1_s = k1_s + (k1_s < 0 ? K1 : 0);
				k2_s = k2-8;
				k2_s = k2_s + (k2_s < 0 ? K2 : 0);
				k3_s = k3-8;
				k3_s = k3_s + (k3_s < 0 ? K3 : 0);
				//printf("\n");
				
				//consider 26 closest neighbours (+2/-2 in each direction) (later must be made cutoff dependent so that the same "distance" is taken regardless of the grid parameters)
				for (i=0; i<17; i++) {
					k1_sn = k1_s+i;
					k1_sn = k1_sn - (k1_sn >= K1 ? K1 : 0);
					for (j=0; j<17; j++) {
						k2_sn = k2_s+j;
						k2_sn = k2_sn - (k2_sn >= K2 ? K2 : 0);
						for (l=0; l<17; l++) {
							k3_sn = k3_s+l;
							k3_sn = k3_sn - (k3_sn >= K3 ? K3 : 0);

							//calculate distance of neighbour
							dist_neighbour= pow(pow(abs(8-i)*cell[3]/(K1-1),2)+pow(abs(8-j)*cell[7]/(K2-1),2)+pow(abs(8-l)*cell[11]/(K3-1),2),0.5);							
							//printf(" -distance neighbour:%f\n", dist_neighbour);
							
							//make sure we ignore the current grid point where the potential is being calculated
							if ( k1_sn!=k1 || k2_sn!=k2 || k3_sn!=k3 ) {
							
								C=100*i_4pi_e0/dist_neighbour;
								//printf(" -C1:%f\n", C);
										
								//calculate neighbour index
								ind_neighbour=k1_sn*dim2*dim3+k2_sn*dim3+k3_sn;
								//printf(" -ind_neighbour:%d\n", ind_neighbour);
										
								//calculate potential due to current neighbour
								C*=q_arr_charge[ind_neighbour]; //*erfc(ewald*dist_neighbour); //for debugging purpose: we remove the erfc, makes it easier to obtain a smooth visualisation
								//printf(" -C2:%1.16f\n", C);
								/*if ( k1==0 && k2==0 && k3==0 ) {
									printf("ind_current: %d, k1_sn: %d, k2_sn: %d, k3_sn: %d, ind_neighbour: %d, distance: %1.16f, potential contribution: %1.16f\n", ind_current, k1_sn, k2_sn, k3_sn, ind_neighbour, dist_neighbour, C);
								}
								if ( k1==0 && k2==0 && k3==1 ) {
									printf("ind_current: %d, k1_sn: %d, k2_sn: %d, k3_sn: %d, ind_neighbour: %d, distance: %1.16f, potential contribution: %1.16f\n", ind_current, k1_sn, k2_sn, k3_sn, ind_neighbour, dist_neighbour, C);
								}*/
										
								//update potential
								q_arr_direct[ind_current]+=C;
								//printf(" -cur potential:%1.16f\n", q_arr_direct[ind_current]);											
	
								//debug
								/*if (k1==0) {
									if (k2==0) {
										if (k3==0) {
											printf("  -index neighbour:%d, charge:%1.16f, distance:%1.16f, pot current:%1.16f,\n", ind_neighbour, q_arr_charge[ind_neighbour], dist_neighbour, C);
										}
										if (k3==1) {
											printf("  -index neighbour:%d, charge:%1.16f, distance:%1.16f, pot current:%1.16f,\n", ind_neighbour, q_arr_charge[ind_neighbour], dist_neighbour, C);
										}
									}
								}*/							
								//if ( k1==0 && k2==0 && k3==0) printf("taking into account i=%d, j= %d, l=%d\n", i, j, l);
							} /*else {
								if ( k1==0 && k2==0 && k3==0 ) printf("ignoring i=%d, j= %d, l=%d\n", i, j, l);
								if ( k1==0 && k2==0 && k3==1 ) printf("ignoring i=%d, j= %d, l=%d\n", i, j, l);
							}*/
						}
					}
				}
			}
				//printf("index:%d, direct:%1.16f\n", ind_current,q_arr_direct[ind_current]);
				//counter++;
				//printf("counter=%d\n", counter);
		}
	}
  }
  
  
  return 0;
}

