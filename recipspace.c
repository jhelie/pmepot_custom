/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/* $Id */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#if !defined(M_PI)
#define M_PI 3.14159265358979323846
#endif

static void dftmod(float *bsp_mod, float *bsp_arr, int nfft) {
  int j, k;
  float twopi, arg, sum1, sum2;
  float infft = 1.0/nfft;
/* Computes the modulus of the discrete fourier transform of bsp_arr, */
/*  storing it into bsp_mod */
  twopi =  2.0 * M_PI;

  for (k = 0; k <nfft; ++k) {
    sum1 = 0.;
    sum2 = 0.;
    for (j = 0; j < nfft; ++j) {
      arg = twopi * k * j * infft;
      sum1 += bsp_arr[j] * cos(arg);
      sum2 += bsp_arr[j] * sin(arg);
    }
    bsp_mod[k] = sum1*sum1 + sum2*sum2;
  }
}

void compute_b_spline(float *frac, float *M);

static void compute_b_moduli(float *bm, int K) {
  int i;
  float fr[3];
  
  float M[12];
  float *scratch = malloc(sizeof(float)*K);

  fr[0]=fr[1]=fr[2]=0.0;
  compute_b_spline(fr,M);  
  for (i=0; i<4; i++) bm[i] = M[i];
  for (i=4; i<K; i++) bm[i] = 0.0;
  dftmod(scratch, bm, K);
  for (i=0; i<K; i++) bm[i] = 1.0/scratch[i];

  free(scratch);
}

//exp(-(pi/alpha)^2*b_vector_dim^2*n^2) for n between 0 and nb_grid_point_dim/2
static void init_exp(float *xp, int K, float recip, float piob) {
  int i;
  float fac;
  fac = -piob*recip*recip;
  //why only K/2?
  //because (apparently): the expression would be between -K/2 and K/2 and because there's a square involved xp[-i]=xp[i]
  //this probably has something to do with the centering of the origin in the reciprocal space (with the -K/2 translation)
  for (i=0; i<= K/2; i++)
    xp[i] = exp(fac*i*i);
} 

static float dot_product(const float *v1, const float *v2) {
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

static void cross_product(float *c, const float *v1, const float *v2) {
  c[0] = v1[1]*v2[2]-v2[1]*v1[2];
  c[1] = v2[0]*v1[2]-v1[0]*v2[2];
  c[2] = v1[0]*v2[1]-v2[0]*v1[1];
}

static float cell_volume(const float *cell) {
  float v, t[3];
  cross_product(t,&cell[3],&cell[6]);
  v = dot_product(t,&cell[9]);
  return abs(v);
}

static int is_parallel(const float *v1, const float *v2) {
  float t[3];
  cross_product(t,v1,v2);
  return ( t[0] == 0. && t[1] == 0. && t[2] == 0. );
}


//This function is called by the pmepot_add function for each frame
//*****************************************************************
float compute_energy_recip(float *q_arr_recip, const float *cell, const float *rcell, const int *dims, float ewald) {

  //NB: this function actually updates the value of q_arr_recip so that they should contain the fourier coefficient of the potential?

  float energy = 0.0;
  float piob, i_pi_volume;
  float *bm1, *bm2, *bm3, *exp1, *exp2, *exp3;
  int pad2, pad3;
  int k1, k2, k3, ind;
  int K1, K2, K3, dim2, dim3;

  //stores number of grid points in each direction in K1 (gridA), K2(gridB) and K3(gridC)
  K1=dims[0]; K2=dims[1]; K3=dims[2]; dim2= dims[3]; dim3=dims[4];
  pad2 = (dim2-K2)*dim3;
  pad3 = dim3-K3-(K3 & 1 ? 1 : 2);

  //list of numbers of incremental modulo in the reciprocal space n1/K1 for n1 between 0 and K1?
  bm1 = malloc(sizeof(float)*K1);
  bm2 = malloc(sizeof(float)*K2);
  bm3 = malloc(sizeof(float)*K3);

  //each exponential is a list of K/2 float numbers: exp(-(pi/alpha)^2*b_vector^2*n^2) for n between 0 and K/2
  exp1 = malloc(sizeof(float)*(K1/2 + 1));
  exp2 = malloc(sizeof(float)*(K2/2 + 1));
  exp3 = malloc(sizeof(float)*(K3/2 + 1));

  //stores 1 / the modulus of the discrete fourier transform of bsp_arr . NB: could this be to nornalise the DFT involved in the charges interpolation calculation?
  compute_b_moduli(bm1, K1);
  compute_b_moduli(bm2, K2);
  compute_b_moduli(bm3, K3);

  //stores the value of the term 1/(Pi*L3) where L is the box dimension
  i_pi_volume = 1.0/(M_PI * cell_volume(cell));
  
  //stores the value of (Pi/alpha)^2
  piob = M_PI/ewald;
  piob *= piob;

  /* for (n=0; n<6; virial[n++] = 0.0); */

  //case where the cell is a cuboid. i.e. A, B and C vectors are all orthogonal to each other
  if ( cell[4] == 0. && cell[5] == 0. && cell[6] == 0. && cell[8] == 0. && cell[9] == 0. && cell[10] == 0. ) {

    //retrieve the non null coordinates of the reciprocal vectors (same as the direct vectors)
    float recipx = rcell[3];
    float recipy = rcell[7];
    float recipz = rcell[11];
    
    //calculate exp(-(pi/alpha)^2*b_vector^2*n^2) for 0<n<K/2 for K1, K2 and K3 
    //K=nb of grid points in each dimension
    //-> this is prob "where" the energy is calculated, i.e. grid point in each dimension
    init_exp(exp1, K1, recipx, piob);
    init_exp(exp2, K2, recipy, piob);
    init_exp(exp3, K3, recipz, piob);

    //sum over each dimensions
    ind = 0;
    for ( k1=0; k1<K1; ++k1 ) {				//for each grid point in the K1 dimension (gridA)
      float m1, m11, b1, xp1;
      int k1_s = k1<=K1/2 ? k1 : k1-K1;		//past K1/2 we consider the corresponding negative value: efficient way to cover -K1/2:K1/2 range while indexing from 0 to K1
      b1 = bm1[k1];							//should it not be k1_s?
      m1 = k1_s*recipx;						//position on the grid
      m11 = m1*m1;							//square of that
      xp1 = i_pi_volume*exp1[abs(k1_s)];	// 1/(Pi*L^3) * exp[- (pi/alpha)^2 * recipx^2 * (k1_s)^2]
      for ( k2=0; k2<K2; ++k2 ) {
        float m2, m22, b1b2, xp2;
        int k2_s = k2<=K2/2 ? k2 : k2-K2;
        b1b2 = b1*bm2[k2];
        m2 = k2_s*recipy;
        m22 = m2*m2;
        xp2 = exp2[abs(k2_s)]*xp1;			// 1/(Pi*L^3) * exp[ -(pi/alpha)^2 (recipx^2 * (k1_s)^2 + recipy^2 * (k2_s)^2)]
        
        //this test excludes the contribution from the origin and avoids a division by 0 when calculating imsq
        if ( k1==0 && k2==0 ) {
          q_arr_recip[ind++] = 0.0;
          q_arr_recip[ind++] = 0.0;
          k3 = 1;
        } else {
          k3 = 0;
        }
        
        for ( ; k3<=K3/2; ++k3 ) {			// why only K/2 for this dimension? => it could be due to the fact that both q_arr_recip[ind] and q_arr_recip[ind+1] are updated at the same time.
          float m3, m33, xp3, msq, imsq, vir, fac;	//i.e. ind corresponds to k3 and ind+1 to -k3?s
          float theta3, theta, q2, qr, qc, C;
          theta3 = bm3[k3] *b1b2;			// inverse of the module of the DFT
          m3 = k3*recipz;					// position on the grid
          m33 = m3*m3;						// square of that
          msq = m11 + m22 + m33;			// sum of the square position on the grid
          imsq = 1.0/msq;					// its inverse
          xp3 = exp3[k3];					// 1/(Pi*L^3) * exp[- (pi/alpha)^2 * recipz^2 * k3^2]
          C = xp2*xp3*imsq;					// 1/imsq * 1/(Pi*L^3) * exp[ -(pi/alpha)^2 * ( recipx^2*(k1_s)^2 + recipy^2*(k2_s)^2 + recipz^2*k3^2 ) ]
          theta = theta3*C;          


          //only useful to calculate the variable energy
          qr = q_arr_recip[ind];
          qc = q_arr_recip[ind+1];
          q2 = 2*(qr*qr + qc*qc)*theta3;
          if ( (k3 == 0) || ( k3 == K3/2 && ! (K3 & 1) ) ) q2 *= 0.5;
		  //

          //*******************************
          q_arr_recip[ind] *= theta;			// why modify the values in q_arr_recip? => because the variable energy is actually useless, it's q_arr_recip that is used to store the potential
          q_arr_recip[ind+1] *= theta;			// now the question is what does the index "ind" correspond to?
		  //printf("-index: %d, coef: %1.16f\n", ind, q_arr_recip[ind]);
          //*******************************
          vir = -2*(piob+imsq);
          fac = q2*C;						//
          energy += fac;					// the energy appears to be a total energy, not an array with a position
          /*
          virial[0] += fac*(1.0+vir*m11);
          virial[1] += fac*vir*m1*m2;
          virial[2] += fac*vir*m1*m3;
          virial[3] += fac*(1.0+vir*m22);
          virial[4] += fac*vir*m2*m3;
          virial[5] += fac*(1.0+vir*m33);
          */
          ind += 2;
        }
        ind += pad3;
      }
      ind += pad2;
    }

  //case where the cell is NOT a cuboid (ignored for now) but still some regularity?
  } else if ( is_parallel(&cell[9],&rcell[9]) ) {

    const float *recip1, *recip2, *recip3;
    float recip3_x, recip3_y, recip3_z;
    recip1 = &rcell[3];
    recip2 = &rcell[6];
    recip3 = &rcell[9];
    recip3_x = recip3[0];
    recip3_y = recip3[1];
    recip3_z = recip3[2];
    init_exp(exp3, K3, sqrt(dot_product(recip3,recip3)), piob);

    ind = 0;
    for ( k1=0; k1<K1; ++k1 ) {
      float b1, m1[3];
      int k1_s = k1<=K1/2 ? k1 : k1-K1;
      b1 = bm1[k1];
      m1[0] = k1_s*recip1[0];
      m1[1] = k1_s*recip1[1];
      m1[2] = k1_s*recip1[2];
      /* xp1 = i_pi_volume*exp1[abs(k1_s)]; */
      for ( k2=0; k2<K2; ++k2 ) {
        float xp2, b1b2, m2_x, m2_y, m2_z;
        int k2_s = k2<=K2/2 ? k2 : k2-K2;
        b1b2 = b1*bm2[k2];
        m2_x = m1[0] + k2_s*recip2[0];
        m2_y = m1[1] + k2_s*recip2[1];
        m2_z = m1[2] + k2_s*recip2[2];
        /* xp2 = exp2[abs(k2_s)]*xp1; */
        xp2 = i_pi_volume*exp(-piob*(m2_x*m2_x+m2_y*m2_y+m2_z*m2_z));
        if ( k1==0 && k2==0 ) {
          q_arr_recip[ind++] = 0.0;
          q_arr_recip[ind++] = 0.0;
          k3 = 1;
        } else {
          k3 = 0;
        }
        for ( ; k3<=K3/2; ++k3 ) {
          float xp3, msq, imsq, vir, fac;
          float theta3, theta, q2, qr, qc, C;
          float m_x, m_y, m_z;
          theta3 = bm3[k3] *b1b2;
          m_x = m2_x + k3*recip3_x;
          m_y = m2_y + k3*recip3_y;
          m_z = m2_z + k3*recip3_z;
          msq = m_x*m_x + m_y*m_y + m_z*m_z;
          xp3 = exp3[k3];
          qr = q_arr_recip[ind]; qc=q_arr_recip[ind+1];
          q2 = 2*(qr*qr + qc*qc)*theta3;
          if ( (k3 == 0) || ( k3 == K3/2 && ! (K3 & 1) ) ) q2 *= 0.5;
          imsq = 1.0/msq;
          C = xp2*xp3*imsq;
          theta = theta3*C;
          q_arr_recip[ind] *= theta;
          q_arr_recip[ind+1] *= theta;
          vir = -2*(piob+imsq);
          fac = q2*C;
          energy += fac;
          /*
          virial[0] += fac*(1.0+vir*m_x*m_x);
          virial[1] += fac*vir*m_x*m_y;
          virial[2] += fac*vir*m_x*m_z;
          virial[3] += fac*(1.0+vir*m_y*m_y);
          virial[4] += fac*vir*m_y*m_z;
          virial[5] += fac*(1.0+vir*m_z*m_z);
          */
          ind += 2;
        }
        ind += pad3;
      }
      ind += pad2;
    }

  //case for other box geometry
  } else {

    const float *recip1, *recip2, *recip3;
    float recip3_x, recip3_y, recip3_z;
    recip1 = &rcell[3];
    recip2 = &rcell[6];
    recip3 = &rcell[9];
    recip3_x = recip3[0];
    recip3_y = recip3[1];
    recip3_z = recip3[2];

    ind = 0;
    for ( k1=0; k1<K1; ++k1 ) {
      float b1, m1[3];
      int k1_s = k1<=K1/2 ? k1 : k1-K1;
      b1 = bm1[k1];
      m1[0] = k1_s*recip1[0];
      m1[1] = k1_s*recip1[1];
      m1[2] = k1_s*recip1[2];
      /* xp1 = i_pi_volume*exp1[abs(k1_s)]; */
      for ( k2=0; k2<K2; ++k2 ) {
        float b1b2, m2_x, m2_y, m2_z;
        int k2_s = k2<=K2/2 ? k2 : k2-K2;
        b1b2 = b1*bm2[k2];
        m2_x = m1[0] + k2_s*recip2[0];
        m2_y = m1[1] + k2_s*recip2[1];
        m2_z = m1[2] + k2_s*recip2[2];
        /* xp2 = exp2[abs(k2_s)]*xp1; */
        if ( k1==0 && k2==0 ) {
          q_arr_recip[ind++] = 0.0;
          q_arr_recip[ind++] = 0.0;
          k3 = 1;
        } else {
          k3 = 0;
        }
        for ( ; k3<=K3/2; ++k3 ) {
          float xp3, msq, imsq, vir, fac;
          float theta3, theta, q2, qr, qc, C;
          float m_x, m_y, m_z;
          theta3 = bm3[k3] *b1b2;
          m_x = m2_x + k3*recip3_x;
          m_y = m2_y + k3*recip3_y;
          m_z = m2_z + k3*recip3_z;
          msq = m_x*m_x + m_y*m_y + m_z*m_z;
          /* xp3 = exp3[k3]; */
          xp3 = i_pi_volume*exp(-piob*msq);
          qr = q_arr_recip[ind]; qc=q_arr_recip[ind+1];
          q2 = 2*(qr*qr + qc*qc)*theta3;
          if ( (k3 == 0) || ( k3 == K3/2 && ! (K3 & 1) ) ) q2 *= 0.5;
          imsq = 1.0/msq;
          C = xp3*imsq;
          theta = theta3*C;
          q_arr_recip[ind] *= theta;
          q_arr_recip[ind+1] *= theta;
          vir = -2*(piob+imsq);
          fac = q2*C;
          energy += fac;
          /*
          virial[0] += fac*(1.0+vir*m_x*m_x);
          virial[1] += fac*vir*m_x*m_y;
          virial[2] += fac*vir*m_x*m_z;
          virial[3] += fac*(1.0+vir*m_y*m_y);
          virial[4] += fac*vir*m_y*m_z;
          virial[5] += fac*(1.0+vir*m_z*m_z);
          */
          ind += 2;
        }
        ind += pad3;
      }
      ind += pad2;
    }

  }

  free(bm1);
  free(bm2);
  free(bm3);
  
  free(exp1);
  free(exp2);
  free(exp3);

  /* for (n=0; n<6; ++n) virial[n] *= 0.5; */
  return 0.5*energy;
}

