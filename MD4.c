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
    double delta_r_k;
    double gravitation;
    double tmp_pos[3];
    double tmp_fli[Ndim];
        for(i=0;i<Nbody;i++)
    {

      for (l = 0; l < Ndim; l++)
      {
        tmp_pos[l] = pos[l][i];
        tmp_fli[l] = f[l][i];
      }

      for(j=i+1;j<Nbody;j++)
      {
        delta_pos[0][k] = tmp_pos[0] - pos[0][j];
        delta_pos[1][k] = tmp_pos[1] - pos[1][j];
        delta_pos[2][k] = tmp_pos[2] - pos[2][j];

        delta_r_k = sqrt(delta_pos[0][k] * delta_pos[0][k] +
          delta_pos[1][k] * delta_pos[1][k] +
          delta_pos[2][k] * delta_pos[2][k]);
        size = radius[i] + radius[j];
        double G_mass_r3 = G * mass[i] * mass[j] / (delta_r_k * delta_r_k * delta_r_k);
        if (delta_r_k >= size)
        {
          for (l = 0; l < Ndim; l++)
          {
            gravitation = delta_pos[l][k] * G_mass_r3;
            tmp_fli[l] -= gravitation;
            f[l][j] += gravitation;
          }
        }
        else
        {
          for (l = 0; l < Ndim; l++)
          {
            gravitation = delta_pos[l][k] * G_mass_r3;
            tmp_fli[l] += gravitation;
            f[l][j] -= gravitation;
          }
          collisions++;
        }
        k = k + 1;
      }
      for (l = 0; l < Ndim; l++)
        f[l][i] = tmp_fli[l];
        }
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




