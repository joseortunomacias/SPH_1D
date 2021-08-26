/* In order to compile this source code, yo need to compile as:

gcc -fopenmp -o SPH_indexx SPH_indexx.c utility.o my.o -lm

and need in the same folder, the objects utility.o and my.o (see how to obtain them in "makefile")
and the files "utility.h" and "my.h"
*/

#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>
#include <libgen.h>
#include <math.h>

#include <omp.h>

#include "utility.h"
#include "my.h"

#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))


typedef struct{
  double x;
  double v;
  double rho;
  double e;
  double h;
  
  double v_coeff;
  double rho_coeff;
  double e_coeff;
  
} Particle;


//----------------------------------------
// PROTOTYPES
//----------------------------------------
double    dW(double, double, double);	 	 
double    set_h(int, int, int, Particle *);
double    set_h_and_index(int, int, int, Particle *, long *);
int       qsort_r_particles(const void *, const void *);
void      get_coeff(int, int, Particle *, long *);






//----------------------------------------
// MAIN
//----------------------------------------
void main()
{
  Particle  *particles; 
  int        i, k;
  double     interval_x;
  int        N_part, N_steps, N_SPH;
  double     t_step;
  double     h;
  FILE      *fp_tmp, *fp_tmp1, *fp_tmp2;
  long     **index;

  N_part   = 1501;
  N_steps  = 200;
  N_SPH    = 100;
  t_step   = 0.0005;








//initial conditions
  particles  = (Particle *) calloc(N_part, sizeof(Particle));
  interval_x = 1.0/(N_part-1.0);
  index      = (long **) calloc(N_part, sizeof(long *));
  for(i=0; i<N_part; i++) {index[i] = (long *) calloc(N_part, sizeof(long));}


  for(i=0; i<N_part; i++)
  {
    particles[i].x    = interval_x*i;
    particles[i].rho  = 1.0/interval_x;
    particles[i].e    = 1.0e-5;
  }
  
  particles[N_part/2].e = 1.0;




//Time stepping loop
  fp_tmp = fopen("densidad.dat", "w");
  for(k=0; k<N_steps; k++)
  {
    //set_h
    #pragma omp parallel for default(none) private(i) shared(N_part, N_SPH, particles, index) 
    for(i=0; i<N_part; i++)
    {
      particles[i].h = set_h_and_index(i, N_part,  N_SPH,  particles, index[i]);
    }

    //getSPHcoefficients
    #pragma omp parallel for default(none) private(i) shared(N_SPH, N_part, particles, index) 
    for(i=0; i<N_part; i++)
    {
      get_coeff(i, N_SPH, particles, index[i]);
    }

    //Perform integration
    #pragma omp parallel for default(none) private(i) shared(N_part, t_step, particles) 
    for(i=0; i<N_part; i++)
    {
      particles[i].x   = particles[i].x   + (particles[i].v)*t_step;
      particles[i].v   = particles[i].v   + (particles[i].v_coeff)*t_step;
      particles[i].rho = particles[i].rho + (particles[i].rho_coeff)*t_step;
      particles[i].e   = particles[i].e   + (particles[i].e_coeff)*t_step;
    } 
  
    //Writte into a file
    for(i=0; i<N_part; i++){fprintf(fp_tmp,"%d	%f \n",i, particles[i].rho);}
  
  }


  fclose(fp_tmp);
  free(particles);

}



//*******************************************************************************
//----------------------------------------
// dW_ij/dx_i
//----------------------------------------

double   dW(double x_i, double x_j, double h)
{
  double r, r_h;
  double sign;


  sign=1.0;
  if ((x_i - x_j)<0.0){sign=-1.0;}
  r   = sign*(x_i - x_j);//absolute value
  r_h = r/h;  // in order to reduce the number of operations (becuase it will appear several more times)

  if(r_h <= 1.0)       {return ( sign*(-2.0/pow3(h)*r + 3.0/2.0*pow2(r_h)/pow2(h)) );}
  else if(r_h <= 2.0)  {return ( -sign*0.5/pow2(h)*pow2(2.0 - r_h) );}
  else                 {return (0.0);}
}




//----------------------------------------
// set_h
//----------------------------------------
double   set_h(int i, int N_part, int N_SPH, Particle *particles)
{
  double *r_particles;
  double  h;
  int     j;

  r_particles = (double *) calloc(N_part, sizeof(double));

  for(j=0; j<N_part; j++)
  {
    r_particles[j] = fabs(particles[i].x - particles[j].x);
  }
  //sort the distances
  qsort(r_particles, N_part, sizeof(double), qsort_r_particles);

  //h will be the distance to the N_SPH th particle
  h = r_particles[N_SPH];
  free(r_particles);
  return (h);  
}





//----------------------------------------
// qsort_r_particles()
//----------------------------------------
int qsort_r_particles(const void *arg1, const void *arg2)
{
  double  *r_part1, *r_part2;

  r_part1 = (double *) arg1;
  r_part2 = (double *) arg2;  

  if((*r_part1) > (*r_part2))       {return(+1);}
  else if ((*r_part1) < (*r_part2)) {return(-1);}
  else                              {return(0);}
}



//----------------------------------------
// set_h_and_index
//----------------------------------------
double    set_h_and_index(int i, int N_part, int N_SPH, Particle *particles, long *index)
{
  int    j;
  double h, *r_particles;
  

  r_particles = (double *) calloc(N_part, sizeof(double));

  for(j=0; j<N_part; j++)
  {
    r_particles[j] = fabs(particles[i].x - particles[j].x);
  }

  indexx((long)N_part,r_particles-1, index-1 );
  h = r_particles[index[N_SPH]-1];
  return(h);
}




//----------------------------------------
// get_coeff_i()
//----------------------------------------

void    get_coeff(int i, int N_SPH, Particle *particles,long *index)
{
  double  v_coeff, rho_coeff, e_coeff;
  double  pressure_i, pressure_j;
  double  dW_ij, dW_ii;
  int     j;

  v_coeff    = 0.0;
  rho_coeff  = 0.0;
  e_coeff    = 0.0;

  pressure_i = (5.0/3.0 - 1.0)*(particles[i].rho)*(particles[i].e);

  for(j=0; j<=N_SPH; j++)
  {
    pressure_j = (5.0/3.0 - 1.0)*(particles[index[j+1]-1].rho)*(particles[index[j+1]-1].e);
    dW_ij      = dW(particles[i].x, particles[index[j+1]-1].x, particles[i].h);

    v_coeff   += (pressure_j/pow2(particles[index[j+1]-1].rho) + pressure_i/pow2(particles[i].rho))*dW_ij;
    rho_coeff += (particles[i].v - particles[index[j+1]-1].v)*dW_ij;
    e_coeff   += (pressure_j/pow2(particles[index[j+1]-1].rho) + pressure_i/pow2(particles[i].rho))*(particles[i].v - particles[index[j+1]-1].v)*dW_ij;
  }

  dW_ii     = dW(particles[i].x, particles[i].x, particles[i].h);
  v_coeff   = v_coeff - (2.0*pressure_i/pow2(particles[i].rho) )*dW_ii;

  particles[i].v_coeff   = -v_coeff;
  particles[i].rho_coeff = rho_coeff;
  particles[i].e_coeff   = e_coeff/2.0;


}






















