/*
 * routines for computing charges on a grid 
 * 
 * $Id: gridcharges.c,v 1.3 2005/07/20 15:37:39 johns Exp $
 *
 */

#include <string.h>
#include <math.h>
#include <stdio.h>

void compute_b_spline(float *frac, float *M) {
  int j;
  float x,y,z,x1,y1,z1, div;
  float *Mx, *My, *Mz;
  Mx=M-1; My=M+4-1; Mz=M+2*4-1;
  x=frac[0];
  y=frac[1];
  z=frac[2];
  x1=1.0-x; y1=1.0-y; z1=1.0-z;
  /* Do n=3 case first */
  Mx[1]=.5*x1*x1;
  Mx[2]=x1*x + .5;
  Mx[3]=0.5*x*x;
  Mx[4]=0.0;
  My[1]=.5*y1*y1;
  My[2]=y1*y + .5;
  My[3]=0.5*y*y;
  My[4]=0.0;
  Mz[1]=.5*z1*z1;
  Mz[2]=z1*z + .5;
  Mz[3]=0.5*z*z;
  Mz[4]=0.0;
  /* Now finish the job!    */
  div=1.0/(4-1);
  Mx[4] = x*div*Mx[4-1];
  My[4] = y*div*My[4-1];
  Mz[4] = z*div*Mz[4-1];
  for (j=1; j<=4-2; j++) {
    Mx[4-j] = ((x+j)*Mx[4-j-1] + (4-x-j)*Mx[4-j])*div;
    My[4-j] = ((y+j)*My[4-j-1] + (4-y-j)*My[4-j])*div;
    Mz[4-j] = ((z+j)*Mz[4-j-1] + (4-z-j)*Mz[4-j])*div;
  }
  Mx[1] *= (1.0-x)*div;
  My[1] *= (1.0-y)*div;
  Mz[1] *= (1.0-z)*div;
}

static float dot_product(const float *v1, const float *v2) {
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

static void scale_vector(float *v, float s) {
  v[0] *= s;
  v[1] *= s;
  v[2] *= s;
}

static void scale_add_vector(float *c, float s, const float *v) {
  c[0] += s * v[0];
  c[1] += s * v[1];
  c[2] += s * v[2];
}

static void copy_vector(float *c, const float *v) {
  c[0] = v[0];
  c[1] = v[1];
  c[2] = v[2];
}

static void cross_product(float *c, const float *v1, const float *v2) {
  c[0] = v1[1]*v2[2]-v2[1]*v1[2];
  c[1] = v2[0]*v1[2]-v1[0]*v2[2];
  c[2] = v1[0]*v2[1]-v2[0]*v1[1];
}

void reciprocal_lattice(const float *cell, float *rcell) {
  const float *a1, *a2, *a3;
  float *b1, *b2, *b3;
  //copy coordinates of origin (note: cell=origin + 3 vectors, i.e. 12 coords)
  rcell[0] = cell[0]; rcell[1] = cell[1]; rcell[2] = cell[2];
  
  //define direct and reciprocal vectors as set of 3 coords
  a1 = &cell[1*3]; a2 = &cell[2*3]; a3 = &cell[3*3];
  b1 = &rcell[1*3]; b2 = &rcell[2*3]; b3 = &rcell[3*3];
        
  //now calculate those reciprocal coords
  cross_product(b1,a2,a3);  scale_vector(b1, 1./dot_product(a1,b1));
  cross_product(b2,a3,a1);  scale_vector(b2, 1./dot_product(a2,b2));
  cross_product(b3,a1,a2);  scale_vector(b3, 1./dot_product(a3,b3));
}

//int fill_charges(const int *dims, const float *cell, int natoms, const float *xyzq, float *q_arr_charge, float *q_arr_recip, float *rcell, float *oddd) {
int fill_charges(const int *dims, const float *cell, int natoms, const float *xyzq, float *q_arr_charge, float *q_arr_recip, float *q_arr_direct, float *rcell, float *oddd) {
  //called by: fill_charges(data->dims,cell,natoms,atoms,q_arr_charge,q_arr_recip,q_arr_direct,rcell,data->oddd);
  
  int i, j, k, l;
  int K1, K2, K3, dim2, dim3;
  float frac[3], Mi[12];
  float ox,oy,oz,r1x,r1y,r1z,r2x,r2y,r2z,r3x,r3y,r3z,kx,ky,kz;

  //stores numbers of grid point in each direction in K1 (gridA), K2 (gridB) and K3 (gridC) 
  K1=dims[0]; K2=dims[1]; K3=dims[2];
  
  //dim2=gridB, dim3=gridC+2
  dim2=dims[3]; dim3=dims[4];

  //initializes the q_arr_recip with 0 values 
  memset( (void*) q_arr_charge, 0, K1*dim2*dim3 * sizeof(float) );
  memset( (void*) q_arr_recip, 0, K1*dim2*dim3 * sizeof(float) );
  memset( (void*) q_arr_direct, 0, K1*dim2*dim3 * sizeof(float) );

  //calculate reciprocal cell vectors (b1,b2,b3) from (a1,a2,a3) (which corresponds to A,B and C vectors in the GUI)
  reciprocal_lattice(cell,rcell);
  
  //center the reciprocal cell origin by removing half a direct cell vector along each coordinate
  //NB: shouldn't it be rcell[0],rcell[1] and rcell[2]? -> after trying it apparently not...
  scale_add_vector(&rcell[0], -0.5*(K1-1.)/K1, &cell[1*3]);
  scale_add_vector(&rcell[0], -0.5*(K2-1.)/K2, &cell[2*3]);
  scale_add_vector(&rcell[0], -0.5*(K3-1.)/K3, &cell[3*3]);

  //copy reciprocal vectors into oddd vectors as well as updated reciprocal origin
  copy_vector(&oddd[0*3],&rcell[0*3]);
  copy_vector(&oddd[1*3],&cell[1*3]);
  copy_vector(&oddd[2*3],&cell[2*3]);
  copy_vector(&oddd[3*3],&cell[3*3]);
  
  //scale the oddd vectors by the number of grid points in each directions
  //this means that 1 reciprocal grid point + 1 reciprocal unit vector = next reciprocal grid point in that direction
  scale_vector(&oddd[1*3],1./K1);
  scale_vector(&oddd[2*3],1./K2);
  scale_vector(&oddd[3*3],1./K3);

  //NB: in case the cell is a cuboid only r1x, r2y and r3z are non null
  ox = rcell[0];
  oy = rcell[1];
  oz = rcell[2];
  r1x = rcell[3];
  r1y = rcell[4];
  r1z = rcell[5];
  r2x = rcell[6];
  r2y = rcell[7];
  r2z = rcell[8];
  r3x = rcell[9];
  r3y = rcell[10];
  r3z = rcell[11];
  kx = 2./K1;  /* cancels position shifts below */
  ky = 2./K2;
  kz = 2./K3;

  for (i=0; i<natoms; i++) {
    float px,py,pz,sx,sy,sz;
    float x,y,z,q;
    int u1, u2, u2i, u3i;

	/*printf("atom nb:%d\n", i);
	printf("direct vectors:\n");
	printf(" -a1x=%f\n", cell[3]);
	printf(" -a1y=%f\n", cell[4]);  
	printf(" -a1z=%f\n", cell[5]);  
	printf(" -a2x=%f\n", cell[6]);
	printf(" -a2y=%f\n", cell[7]);  
	printf(" -a2z=%f\n", cell[8]);  
	printf(" -a3x=%f\n", cell[9]);
	printf(" -a3y=%f\n", cell[10]);  
	printf(" -a3z=%f\n", cell[11]);
	printf("reciprocal vectors:\n");
	printf(" -b1x=%f\n", rcell[3]);
	printf(" -b1y=%f\n", rcell[4]);  
	printf(" -b1z=%f\n", rcell[5]);  
	printf(" -b2x=%f\n", rcell[6]);
	printf(" -b2y=%f\n", rcell[7]);  
	printf(" -b2z=%f\n", rcell[8]);  
	printf(" -b3x=%f\n", rcell[9]);
	printf(" -b3y=%f\n", rcell[10]);  
	printf(" -b3z=%f\n", rcell[11]);*/

    //get coordinates of current atom and center them on the origin
    px = xyzq[4*i+0] - ox;
    py = xyzq[4*i+1] - oy;
    pz = xyzq[4*i+2] - oz;
    
	/*printf(" -px:%f\n", px);
	printf(" -py:%f\n", py);
    printf(" -pz:%f\n", pz);*/

    sx = px*r1x + py*r1y + pz*r1z + kx; 	//in practice = px*r1x+kx
    sy = px*r2x + py*r2y + pz*r2z + ky;		//in practice = py*r2y+ky
    sz = px*r3x + py*r3y + pz*r3z + kz;		//in practice = pz*r3z+kz
    
    /*printf(" -kx=%f\n", kx);
    printf(" -ky=%f\n", ky);
    printf(" -kz=%f\n", kz);
    printf(" -r1x=%f\n", r1x);
    printf(" -r1y=%f\n", r1y);
    printf(" -r1z=%f\n", r1z);
    printf(" -r2x=%f\n", r2x);    
    printf(" -r2y=%f\n", r2y);
    printf(" -r2z=%f\n", r2z);
    printf(" -r3x=%f\n", r3x);
    printf(" -r3y=%f\n", r3y);   
    printf(" -r3z=%f\n", r3z);
    printf(" -sx=%f\n", sx);
    printf(" -sy=%f\n", sy);
    printf(" -sz=%f\n", sz);*/
    
    x = K1 * ( sx - floor(sx) );
    y = K2 * ( sy - floor(sy) );
    z = K3 * ( sz - floor(sz) );
    
    /*printf(" -x=%f\n", x);
    printf(" -y=%f\n", y);
    printf(" -z=%f\n", z);*/
    /*  Check for rare rounding condition where K * ( 1 - epsilon ) == K */
    /*  which was observed with g++ on Intel x86 architecture.           */
    if ( x == K1 ) x = 0;
    if ( y == K2 ) y = 0;
    if ( z == K3 ) z = 0;

    q = xyzq[4*i+3];

	//printf(" -charge:%f\n", q);


    //these lines determine the closest point on the grid (u1, u2i, u3i) corresponding to the atom coordinates (with the convention that the lowest index node is chosen)
    //-------------------------------------------------------------------------------------------------------------------------------------------------------------------
    u1 = (int)x;		//si l'on se rappelle que Ki correspond au nombre de point dans une direction
    u2i = (int)y;		//ui correspond au point avant la ou est situe x
    u3i = (int)z;		//sachant que x evolue entre 0 et K1, les positions x sont symetrises par rapport a K/2
    /*printf(" -u1=%d\n", u1);
    printf(" -u2i=%d\n", u2i);
    printf(" -u3i=%d\n", u3i);    */
    frac[0] = x - u1;
    frac[1] = y - u2i;
    frac[2] = z - u3i;
    compute_b_spline(frac,Mi);
    u1 -= 4;
    u2i -= 4;
    u3i -= 4;
    u3i++;
	//not sure exactly how grid points are accessed/index from the atoms x,y,z coordinates 
    //but essentially it seems that q_arr_recip follows the dx file convention, i.e. z varies first, then y then x.
    
    //this loop now takes into account the 48 closest neighbours? (3 loops of 4?) 
    //-----------------------------------------------------------
    //NB: the rounding above is just a clever way to calculate the upper values of u1,u2,u3: 4 values are then taken into account ui-3 to u1 and the
    //actual position of the current atom is between ui-2 and ui-1
    //SO IN SUMMARY: for each point, for each dimension the 2 nearest grid points before and after the real position are taken into account
    for (j=0; j<4; j++) {
      //printf("  -j=%d\n", j);
      float m1;
      int ind1;
      m1 = Mi[j]*q;			//scale the charge of the current atom based on x deviation to the current grid point
      u1++;
      ind1 = u1 + (u1 < 0 ? K1 : 0); //this test takes into account periodic boundary conditions!
      //printf("  -ind1=%d\n", ind1);
      u2 = u2i;
      for (k=0; k<4; k++) {
        //printf("   -k=%d\n", k);
        float m1m2;
		int ind2;
        m1m2 = m1*Mi[4+k];	//scale it again based on y deviation to the current grid point
		u2++;
		ind2 = ind1*dim2 + (u2 + (u2 < 0 ? K2 : 0));
        //printf("   -ind2=%d\n", ind2);
        for (l=0; l<4; l++) {
		  //printf("    -l=%d\n", l);
		  float m3;
		  int ind;
		  int u3 = u3i + l;
		  //printf("    -u3=%d\n", u3);
		  m3 = Mi[2*4 + l]; //and finally scale it based on deviation to the current grid point
          ind = ind2*dim3 + (u3 + (u3 < 0 ? K3 : 0));
          //printf("     -ind=%d\n", ind);
          
          q_arr_charge[ind] += m1m2*m3;	//if q_arr_charge and q_arr_recip point to the same variable we're effectively multiplying it by two
          q_arr_recip[ind] += m1m2*m3;
		  q_arr_direct[ind] = 0.0;
        }
      }
    }
  }

  return 0;
}

