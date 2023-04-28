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
      for(int i=0;i<Nbody;i++)
        r[i] = 1.0/mass[i];
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
    double tmp_pos0;
    double tmp_pos1;
    double tmp_pos2;
    double tmp_fli0;
    double tmp_fli1;
    double tmp_fli2;
    for(i=0;i<Nbody;i++)
    {

tmp_pos0 = pos[0][i];
tmp_pos1 = pos[1][i];
tmp_pos2 = pos[2][i];


tmp_fli0 = f[0][i];
tmp_fli1 = f[1][i];
tmp_fli2 = f[2][i];
      for(j=i+1;j<Nbody;j++)
      {
        delta_pos[0][k] = tmp_pos0 - pos[0][j];
        delta_pos[1][k] = tmp_pos1 - pos[1][j];
        delta_pos[2][k] = tmp_pos2 - pos[2][j];

        delta_r_k = sqrt(delta_pos[0][k] * delta_pos[0][k] +
          delta_pos[1][k] * delta_pos[1][k] +
          delta_pos[2][k] * delta_pos[2][k]);
        size = radius[i] + radius[j];
        double G_mass_r3 = G * mass[i] * mass[j] / (delta_r_k * delta_r_k * delta_r_k);

            tmp_fli0 -= delta_pos[0][k] * G_mass_r3;
            f[0][j] += delta_pos[0][k] * G_mass_r3;

            tmp_fli1 -= delta_pos[1][k] * G_mass_r3;
            f[1][j] += delta_pos[1][k] * G_mass_r3;


            tmp_fli2-= delta_pos[2][k] * G_mass_r3;
            f[2][j] += delta_pos[2][k] * G_mass_r3;
#if 0
        if (delta_r_k >= size)
        {
            tmp_fli[0] -= delta_pos[0][k] * G_mass_r3;
            f[0][j] += delta_pos[0][k] * G_mass_r3;

            tmp_fli[1] -= delta_pos[1][k] * G_mass_r3;
            f[1][j] += delta_pos[1][k] * G_mass_r3;


            tmp_fli[2] -= delta_pos[2][k] * G_mass_r3;
            f[2][j] += delta_pos[2][k] * G_mass_r3;
          //for (l = 0; l < Ndim; l++)
          //{
            //gravitation = delta_pos[l][k] * G_mass_r3;
            //tmp_fli[l] -= delta_pos[l][k] * G_mass_r3;
            //f[l][j] += delta_pos[l][k] * G_mass_r3;
          //}
        }
        else
        {
            tmp_fli[0] += delta_pos[0][k] * G_mass_r3;
            f[0][j] -= delta_pos[0][k] * G_mass_r3;

            tmp_fli[1] += delta_pos[1][k] * G_mass_r3;
            f[1][j] -= delta_pos[1][k] * G_mass_r3;


            tmp_fli[2] += delta_pos[2][k] * G_mass_r3;
            f[2][j] -= delta_pos[2][k] * G_mass_r3;
          collisions++;
        }
#endif
        k = k + 1;
      }
      f[0][i] = tmp_fli0;
      f[1][i] = tmp_fli1;
      f[2][i] = tmp_fli2;
    }
    for(i=0;i<Nbody;i++)
    {
      pos[0][i] =  pos [0][i] + dt * velo[0][i];
      velo[0][i] = velo[0][i] + dt * (f[0][i] *r[i]);
      
      pos[1][i] = pos[1][i] + dt * velo[1][i];
      velo[1][i] = velo[1][i] + dt * (f[1][i]*r[i]);

      pos[2][i] = pos[2][i] + dt * velo[2][i];
      velo[2][i] = velo[2][i] + dt * (f[2][i] *r[i]);
    }


      }

}




