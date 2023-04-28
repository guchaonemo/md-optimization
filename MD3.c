/*
 *  Simple molecular dynamics code.
 *  2022
 */
#include <stdio.h>
#include <math.h>
#include "coord.h"

void vis_forces(int N,double *f, double *vis, double *vel);
void add_norms(int N,double *r, double *delta);
double forces(double W, double delta, double r);
void wind_forces(int N,double *f, double *vis, double vel);

void evolve(int count,double dt){
int step;
int i,j,k,l;
int have_collided;
double size;
/*
 * Loop over timesteps.
 */
      for(step = 1;step<=count;step++){
        printf("timestep %d\n",step);
        printf("collisions %d\n",collisions);

    /* calculate central force */
    for (i = 0; i < Nbody; i++)
    {
        double tmp = sqrt(pos[0][i] * pos[0][i] + pos[1][i] * pos[1][i] + pos[2][i] * pos[2][i]);
        tmp = G * mass[i] * M_central / (tmp*tmp*tmp);
        f[0][i] = -vis[i] * (velo[0][i] + wind[0]) - tmp* pos[0][i];
        f[1][i] = -vis[i] * (velo[1][i] + wind[1]) - tmp* pos[1][i];
        f[2][i] = -vis[i] * (velo[2][i] + wind[2]) - tmp* pos[2][i];
    }
/* calculate pairwise separation of the particles */
        k = 0;
        for(i=0;i<Nbody;i++){
          for(j=i+1;j<Nbody;j++){
            for(l=0;l<Ndim;l++){
              delta_pos[l][k] = pos[l][i] - pos[l][j];
            }
            k = k + 1;
          }
        }

/* calculate norm of separation vector */
        for(k=0;k<Npair;k++){
          delta_r[k] = 0.0;
        }
        for(i=0;i<Ndim;i++){
	  add_norms(Npair,delta_r,delta_pos[i]);
        }
        for(k=0;k<Npair;k++){
          delta_r[k] = sqrt(delta_r[k]);
        }

/*
 * add pairwise forces.
 */
        k = 0;
double gravitation;
for(i=0;i<Nbody;i++)
{
for(j=i+1;j<Nbody;j++)
{
size = radius[i] + radius[j];
double G_mass_r3 = G * mass[i] * mass[j] / (delta_r[k] *
delta_r[k] * delta_r[k]);
if (delta_r[k] >= size)
{
for (l = 0; l < Ndim; l++)
{
gravitation = delta_pos[l][k] * G_mass_r3;
f[l][i] -= gravitation;
f[l][j] += gravitation;
}
}
else
{
for (l = 0; l < Ndim; l++)
{
gravitation = delta_pos[l][k] * G_mass_r3;
f[l][i] += gravitation;
f[l][j] -= gravitation;
}
collisions++;}
k = k + 1;}}

/* update positions */
        for(i=0;i<Nbody;i++){
          for(j=0;j<Ndim;j++){
            pos[j][i] = pos[j][i] + dt * velo[j][i];
          }
        }

/* update velocities */
        for(i=0;i<Nbody;i++){
          for(j=0;j<Ndim;j++){
            velo[j][i] = velo[j][i] + dt * (f[j][i]/mass[i]);
          }
        }


      }

}




