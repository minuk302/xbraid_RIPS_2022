/*BHEADER**********************************************************************
 * Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
 * Produced at the Lawrence Livermore National Laboratory. Written by 
 * Jacob Schroder, Rob Falgout, Tzanio Kolev, Ulrike Yang, Veselin 
 * Dobrev, et al. LLNL-CODE-660355. All rights reserved.
 * 
 * This file is part of XBraid. For support, post issues to the XBraid Github page.
 * 
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (as published by the Free Software
 * Foundation) version 2.1 dated February 1999.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
 * License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc., 59
 * Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 ***********************************************************************EHEADER*/

 /**
 * Example:       trischur-ex-04.c
 *
 * Interface:     C
 * `
 * Requires:      only C-language support     
 *
 * Compile with:  make trischur-ex-04
 *
 * Description:  Solves a simple optimal control problem in time-parallel:
 * 
 *                 min   \int_0^1 u_1(t)^2 + u_2(t)^2 + gamma c(t)^2  dt
 * 
 *                  s.t.  d/dt u_1(t) = u_2(t)
 *                        d/dt u_2(t) = -u_2(t) + c(t)
 * 
 *               with initial condition u_1(0) = 0, u_2(0) = -1
 *               and piecewise constant control c(t).  
 *
 *               Implements a steepest-descent optimization iteration
 *               using fixed step size for design updates.   
 **/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "braid.h"
#include "braid_test.h"

/*--------------------------------------------------------------------------
 * My App and Vector structures
 *--------------------------------------------------------------------------*/

typedef struct _braid_App_struct
{
   int      myid;        /* Rank of the processor */
   int      ntime;       /* Total number of time-steps (starting at time 0) */
   int      ilower;      /* Lower index for my proc */
   int      iupper;      /* Upper index for my proc */
   int      npoints;     /* Number of time points on my proc */
   double **u;           /* Adjoing vectors at each time point on my proc */
   double *input_u;      /* current u  */
   double *input_v;     /* current v */
   double *rhs;         /* rhs */
} my_App;

double u_in = 0;
double u_out = 0;
double tstart, tstop;

/* Define the state vector at one time-step */
typedef struct _braid_Vector_struct
{
   double *values;     /* Holds the R^2 state vector (u_1, u_2) */

} my_Vector;

/*--------------------------------------------------------------------------
 * Vector utility routines
 *--------------------------------------------------------------------------*/

void
vec_create(int size, double **vec_ptr)
{
   *vec_ptr = (double*) malloc( size*sizeof(double) );
}

void
vec_destroy(double *vec)
{
   free(vec);
}

/*------------------------------------*/

void
vec_copy(int size, double *invec, double *outvec)
{
   int i;
   for (i = 0; i < size; i++)
   {
      outvec[i] = invec[i];
   }
}

/*------------------------------------*/

void
vec_axpy(int size, double alpha, double *x, double *y)
{
   int i;
   for (i = 0; i < size; i++)
   {
      y[i] = y[i] + alpha*x[i];
   }
}

/*------------------------------------*/

void
vec_scale(int size, double alpha, double *x)
{
   int i;
   for (i = 0; i < size; i++)
   {
      x[i] = alpha*x[i];
   }
}

void extract(double *r, double *k, double *h, double *g, int v_size) 
{
   for(int i = 0; i < v_size - 1; i++)
   {
      k[i] = r[i];
   }
   for(int i =0; i < v_size; i++)
   {
      h[i] = r[i+v_size - 1];
   }
   for(int i = 0; i < v_size; i++)
   {
      g[i] = r[i + 2 * v_size - 1];
   }
}

void merge(double *r, double *k, double *h, double *g, int v_size){
   for(int i = 0; i < v_size - 1; i++)
   {
      r[i] = k[i];
   }
   for(int i =0; i < v_size; i++)
   {
      r[i+v_size - 1] = h[i];
   }
   for(int i = 0; i < v_size; i++)
   {
      r[i + 2 * v_size - 1] = g[i];
   }
}

/*--------------------------------------------------------------------------
 * KKT component routines
 *--------------------------------------------------------------------------*/

void
apply_Phi(double dt, double *u)
{
}

/*------------------------------------*/

void
apply_PhiAdjoint(double dt, double *w)
{
}

/*------------------------------------*/
//apply_U(dt, utmp, u, v, index);
void
apply_U(double dt, double *jm, double *u, double * v, int index)
{
   //jm[0] = dt * (v[index] * v[index] + v[index+1] * v[index+1])/(u[index] *u[index] * u[index] );
   //if (index < v_size - 1) 
   if (index < 19) 
   {
      jm[0] = dt * (v[index] * v[index] + v[index+1] * v[index+1])/(u[index] *u[index] * u[index] );
   }
   else
   {
      jm[0] = dt * (v[index] * v[index] + v[index+1] * v[index+1])/(u_out * u_out * u_out);
   }
}

/*------------------------------------*/
/* apply_V(dt, INSERT MATRIX HERE, u, v, INSERT INDEX HERE);*/
void
apply_V(double dt, double * rob, double *u, double *v, int index, int v_size, double u_in, double u_out)
{
   if (index == 0)
   {
      rob[0] *=  dt * (1/ u_in + 1/u[0]);
   }
   if (index == v_size - 1)
   {
      rob[0] *=  dt * (1/ u_out + 1/u[index-1]);
   }
   if (index != 0 && index != v_size - 1)
   {
      rob[0] *=  dt * (1/u[index-1]+1/u[index]);
   }
}

/*------------------------------------*/

/* Maps v to w (which may be the same pointer) */
void
apply_D(double dt, double *v, double *w)
{
   w[0] = -v[0]*dt; 
}

void
apply_DAdjoint(double dt, double *v, double *w)
{
   w[0] = -v[0]*dt; 
}

void
apply_Dinv(double dt, double *v, double *w)
{
   w[0] = -v[0]/dt; 
}

/*------------------------------------*/

/* Maps w to v (which may be the same pointer) */
void
apply_DAdjointinv(double dt, double *w, double *v)
{
   double wtmpone = w[0];
   v[0] = wtmpone/dt; 
}


/* -------------------------------------------
 * TriMGRIT wrapper routines
 *--------------------------------------------------------------------------*/

/* Compute A(u) - f */

int my_TriResidual(braid_App       app,
               braid_Vector    uleft,
               braid_Vector    uright,
               braid_Vector    f,
               braid_Vector    r,
               braid_TriStatus status)
{
   double  t, tprev, tnext, dt;
   double *rtmp, *utmp, *htmp, *h_othertemp, *gtmp, *ktmp, *g_othertemp, *u, *v;
   double  *k, *h, *g;
   int     level, index;
   int v_size = app->ntime;
   u = app->input_u;
   v = app->input_v;
   
   
   braid_TriStatusGetTriT(status, &t, &tprev, &tnext);
   braid_TriStatusGetLevel(status, &level);
   braid_TriStatusGetTIndex(status, &index);
   
   /* extract k, h, g from res */
   vec_create(v_size -1, &k);
   vec_create(v_size , &h);
   vec_create(v_size , &g);
   extract(r->values, k, h, g, v_size);

   /* Get the time-step size */
   if (t < tnext)
   {
      dt = tnext - t;
   }
   else
   {
      dt = t - tprev;
   }

   /* Create temporary vectors */
   vec_create(1, &rtmp);
   vec_create(1, &utmp);
   vec_create(1, &h_othertemp);
   vec_create(1, &gtmp);
   vec_create(1, &htmp);
   vec_create(1, &ktmp);
   vec_create(1, &g_othertemp);

   /* rtmp = U_i u */
   vec_copy(1, (r->values), utmp);
   //apply_U(dt, utmp);
   apply_U(dt, utmp, u, v, index);
   vec_copy(1, utmp, rtmp);

   /* rtmp = rtmp + D_i^-T V_i D_i^-1 u */
   vec_copy(1, (r->values), utmp);
   printf("\nthis??\n %d %f", index, utmp[0]);
   apply_Dinv(dt, utmp, utmp);
   apply_V(dt, utmp, u, v, index, v_size, u_in, u_out);
   //apply_V(dt, gamma_1, gamma_2, utmp);
   apply_DAdjointinv( dt, utmp, utmp);
   vec_axpy(1, 1.0, utmp, rtmp);
   /* rtmp = rtmp + Phi_{i+1}^T D_{i+1}^{-T} V_{i+1} D_{i+1}^{-1} Phi_{i+1} u */
   /* This term is zero at time 0, since Phi_0 = 0 */
   if (uright != NULL)
   {
      vec_copy(1, (r->values), utmp);
      apply_Phi(dt, utmp);
      apply_Dinv(dt, utmp, utmp);
      apply_V(dt, utmp, u, v, index, v_size, u_in, u_out);
      //apply_V(dt, gamma_1, gamma_2, utmp);
      apply_DAdjointinv(dt, utmp, utmp);
      apply_PhiAdjoint(dt, utmp);
      vec_axpy(1, 1.0, utmp, rtmp);
   }
   printf("\nthis1??\n %d %f", index, utmp[0]);
   /* Compute action of west block */
   if (uleft != NULL)
   {
      /* rtmp = rtmp - D_i^{-T} V_i D_i^{-1} Phi_i uleft */
      vec_copy(1, (uleft->values), utmp);
      printf("\nthis??\n %d %f", index, utmp[0]);
      apply_Phi(dt, utmp);
      apply_Dinv(dt, utmp, utmp);
      apply_V(dt, utmp, u, v, index, v_size, u_in, u_out);
      //apply_V(dt, gamma_1, gamma_2, utmp);
      apply_DAdjointinv(dt, utmp, utmp);
      vec_axpy(1, -1.0, utmp, rtmp);
   }
   printf("\nthis2??\n %d %f", index, utmp[0]);
   /* Compute action of east block */
   if (uright != NULL)
   {
      /* rtmp = rtmp - Phi_{i+1}^T D_{i+1}^{-T} V_{i+1} D_{i+1}^{-1} uright */
      vec_copy(1, (uright->values), utmp);
      printf("\nthis4??\n %d %f", index, utmp[0]);
      apply_Dinv(dt, utmp, utmp);
      apply_V(dt, utmp, u, v, index, v_size, u_in, u_out);
      //apply_V(dt, gamma_1, gamma_2, utmp);
      apply_DAdjointinv(dt, utmp, utmp);
      apply_PhiAdjoint(dt, utmp);
      vec_axpy(1, -1.0, utmp, rtmp);
   }
   printf("\nthis3??\n %d %f", index, utmp[0]);
   /* subtract rhs kbar (add -k - L^T D^{-T} V D^{-1} g - L^T D^{-T} h */
   htmp[0] = h[index];
   h_othertemp[0] = h[index+1];
   gtmp[0] = g[index];
   g_othertemp[0] = h[index+1];
   ktmp[0] = k[index];
   
   apply_DAdjointinv(dt, htmp, htmp); 
   apply_DAdjointinv(dt, h_othertemp, h_othertemp);
   vec_copy(1,h_othertemp,rtmp);
   vec_axpy(1, -1, htmp, rtmp); 
   
   
   apply_Dinv(dt, gtmp, gtmp);
   apply_V(dt, gtmp, u, v, index, v_size, u_in, u_out);
   apply_DAdjointinv(dt, gtmp, gtmp);
   apply_Dinv(dt, g_othertemp, g_othertemp);
   apply_V(dt, g_othertemp, u, v, index, v_size, u_in, u_out);
   apply_DAdjointinv(dt, g_othertemp, g_othertemp);
   vec_copy(1,g_othertemp,rtmp);
   vec_axpy(1, -1, gtmp, rtmp);

   vec_axpy(1, -1, ktmp, rtmp);
   
   
   if (f!= NULL)
   {
   vec_copy(1, rtmp, (r->values));
   }
   /* Destroy temporary vectors */
   printf("\n\n is it triresidual??? \n %d %f %f", index, rtmp[0], utmp[0]);
   vec_destroy(rtmp);
   vec_destroy(utmp);
   
   return 0;
}   

/*------------------------------------*/



/* approximately solve Solve A(u) = f */
int
my_TriSolve(braid_App       app,
            braid_Vector    uleft,
            braid_Vector    uright,
            braid_Vector    fleft,
            braid_Vector    fright,
            braid_Vector    f,
            braid_Vector    u,
            braid_TriStatus status)
{
   double  t, tprev, tnext, dt;
   double *utmp, *rtmp;
   int index;
   
   /* Get the time-step size */
   braid_TriStatusGetTriT(status, &t, &tprev, &tnext);
   if (t < tnext)
   {
      dt = tnext - t;
   }
   else
   {
      dt = t - tprev;
   }

   /* Create temporary vector */
   vec_create(1, &utmp);

   /* Initialize temporary solution vector */
   vec_copy(1, (u->values), utmp);
   
   /* Compute residual */
   my_TriResidual(app, uleft, uright, f, u, status);

   /* Apply center block preconditioner (multiply by \tilde{C}^-1) to -r
    *
    * Using \tilde{C} = | v_{i}/dt^2 + v_{i+1}/dt^2 + u_{i}|
    *          
    */
   rtmp = (u->values);
   braid_TriStatusGetTIndex(status, &index);
   double * input_u = app->input_u;
   double * input_v = app->input_v;
   int v_size = app->ntime;
   //if (uright != NULL)
   //{
   if (index > 0 && index < v_size - 1)
   {
      rtmp[0] = 1 / (dt) * (1/input_u[index - 1] + 1 / input_u[index]) 
      + 1 / (dt) *(1 / input_u[index] + 1 / input_u[index + 1])
      + dt * ((input_v[index] * input_v[index]) + (input_v[index+1] * input_v[index+1]))/(input_u[index]*input_u[index] * input_u[index]);
   }
   //rtmp[0] = input_v[index]/(dt * dt) + input_v[index+1]/(dt * dt) + input_u[index];
   
   if (index == 0)
   {
      //rtmp[0] = 1 / (dt) *(1/app->u_in + 1/input_u[index]) +1/(dt) *(1/input_u[index] + 1/input_u[index + 1]);
      rtmp[0] = 1 / (dt) * (1/u_in + 1 / u_in) 
      + 1 / (dt) *(1 / input_u[index] + 1 / input_u[index + 1])
      + dt * ((input_v[index] * input_v[index]) + (input_v[index+1] * input_v[index+1]))/(input_u[index]*input_u[index] * input_u[index]);
   }
      
      /*if (index == v_size - 1)
      {
         rtmp[0] = 1/(dt) *(1/input_u[index - 1] + 1/input_u[index]) +1/(dt) *(1/input_u[index] + 1/app->u_out);
      }*/
   /*}
   else
   {
      /* At the rightmost point, use a different center coefficient approximation */
      /*rtmp[0] = input_v[index]/(dt * dt) + input_u[index];
      //rtmp[0] = input_v[index]/(dt * dt) + input_u[index];
      rtmp[0] = 1 / (dt) * (1/input_u[index - 1] + 1 / input_u[index]) 
         + dt * ((input_v[index] * input_v[index]) + (input_v[index+1] * input_v[index+1]))/(input_u[index]*input_u[index] * input_u[index]);
   }*/

   /* Complete residual update */
   vec_axpy(1, 1.0, utmp, (u->values));
   
   /* no refinement */
   braid_TriStatusSetRFactor(status, 1);

   /* Destroy temporary vectors */
   vec_destroy(utmp);

   return 0;
}   

/*------------------------------------*/

/* This is only called from level 0 */

int
my_Init(braid_App     app,
        double        t,
        braid_Vector *u_ptr)
{
   my_Vector *u;

   /* Allocate the vector */
   u = (my_Vector *) malloc(sizeof(my_Vector));
   vec_create(1, &(u->values));

   //u->values[0] = ((double)braid_Rand())/braid_RAND_MAX;
   //u->values[1] = ((double)braid_Rand())/braid_RAND_MAX;
   u->values[0] = 1; // just random value
   *u_ptr = u;
   return 0;
}

/*------------------------------------*/

int
my_Clone(braid_App     app,
         braid_Vector  u,
         braid_Vector *v_ptr)
{
   my_Vector *v;

   /* Allocate the vector */
   v = (my_Vector *) malloc(sizeof(my_Vector));
   vec_create(1, &(v->values));

   /* Clone the values */
   v->values[0] = u->values[0];

   *v_ptr = v;

   return 0;
}

/*------------------------------------*/

int
my_Free(braid_App    app,
        braid_Vector u)
{
   free(u->values);
   free(u);

   return 0;
}

/*------------------------------------*/

int
my_Sum(braid_App     app,
       double        alpha,
       braid_Vector  x,
       double        beta,
       braid_Vector  y)
{

   (y->values)[0] = alpha*(x->values)[0] + beta*(y->values)[0];

   return 0;
}

/*------------------------------------*/

int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
   double dot = (u->values)[0]*(u->values)[0];
   *norm_ptr = sqrt(dot);

   return 0;
}

/*------------------------------------*/

int
my_Access(braid_App          app,
          braid_Vector       u,
          braid_AccessStatus astatus)
{
   int   done, index, ii;

   /* Print solution to file if simulation is over */
   braid_AccessStatusGetDone(astatus, &done);

   if (done)
   {
      braid_AccessStatusGetILowerUpper(astatus, &(app->ilower), &(app->iupper));
      (app->npoints) = (app->iupper) - (app->ilower) + 1;

      /* Allocate w array in app */
      if ((app->u) == NULL)
      {
         (app->u) = (double **) calloc((app->npoints), sizeof(double *));
      }

      braid_AccessStatusGetTIndex(astatus, &index);
      ii = index - (app->ilower);
      if (app->u[ii] != NULL)
      {
         free(app->u[ii]);
      }
      vec_create(1, &(app->u[ii]));
      vec_copy(1, (u->values), (app->u[ii]));
   }
   return 0;
}

/*------------------------------------*/

int
my_BufSize(braid_App           app,
           int                 *size_ptr,
           braid_BufferStatus  bstatus)
{
   *size_ptr = 1*sizeof(double);
   return 0;
}

/*------------------------------------*/

int
my_BufPack(braid_App           app,
           braid_Vector        u,
           void               *buffer,
           braid_BufferStatus  bstatus)
{
   double *dbuffer = buffer;
   int i;

   for(i = 0; i < 1; i++)
   {
      dbuffer[i] = (u->values)[i];
   }

   braid_BufferStatusSetSize( bstatus,  sizeof(double));
   return 0;
}

/*------------------------------------*/

int
my_BufUnpack(braid_App           app,
             void               *buffer,
             braid_Vector       *u_ptr,
             braid_BufferStatus  bstatus)
{
   my_Vector *u = NULL;
   double    *dbuffer = buffer;
   int i;

   /* Allocate memory */
   u = (my_Vector *) malloc(sizeof(my_Vector));
   u->values = (double*) malloc( sizeof(double) );

   /* Unpack the buffer */
   for(i = 0; i < 1; i++)
   {
      (u->values)[i] = dbuffer[i];
   }

   *u_ptr = u;
   return 0;
}


/*solves for e in H_p e = r*/
void solveHp(double *e, double *r, double *u, double *v, double dt, int v_size, double u_in, double u_out){
   braid_Core  core;
   my_App     *app;
      
   int         rank, ntime;
   
   int         max_levels, min_coarse, nrelax, nrelaxc, cfactor, maxiter;
   int         access_level, print_level;
   double tol;
   /* rank should not matter */
   rank = 2;
   /* Define time domain */
   ntime  = v_size;              /* Total number of time-steps */

   /* Define some Braid parameters */
   max_levels     = 30;
   min_coarse     = 1;
   nrelax         = 1;
   nrelaxc        = 7;
   maxiter        = 20;
   cfactor        = 2;
   tol            = 1.0e-6;
   access_level   = 1;
   print_level    = 2;
   /* Set up twohe app structure */
   app = (my_App *) malloc(sizeof(my_App));
   app->myid     = rank;
   app->ntime    = ntime;
   app->u        = NULL;
   app->input_u = u;
   app->input_v = v;
   app->rhs = r;
   /* Initialize k, h, g vectors */
   double * k;
   double * h;
   double * g;
   vec_create(v_size-1, &k);
   vec_create(v_size, &h);
   vec_create(v_size, &g);
   extract(r, k, h, g, v_size);
   /* Initialize XBr-1aid */
   braid_InitTriMGRIT(MPI_COMM_WORLD, MPI_COMM_WORLD, dt, tstop, ntime-1, app,
                      my_TriResidual, my_TriSolve, my_Init, my_Clone, my_Free,
                      my_Sum, my_SpatialNorm, my_Access,
                      my_BufSize, my_BufPack, my_BufUnpack, &core);
 
   /* Set some XBraid(_Adjoint) parameters */
   braid_SetMaxLevels(core, max_levels);
   braid_SetMinCoarse(core, min_coarse);
   braid_SetNRelax(core, -1, nrelax);
  
   if (max_levels > 1)
   {
      braid_SetNRelax(core, max_levels-1, nrelaxc); /* nrelax on coarsest level */
   }
   braid_SetCFactor(core, -1, cfactor);
   braid_SetAccessLevel(core, access_level);

   braid_SetPrintLevel( core, print_level);       
   braid_SetMaxIter(core, maxiter);
   braid_SetAbsTol(core, tol);

   /* Parallel-in-time TriMGRIT simulation */
   braid_Drive(core);

   /*   Compute control v and adjoint w from x */
   double * out_u;
   double * out_v, * out_vtmp, * out_w, *w, * wtmp, *gtmp, *htmp;
   vec_create(ntime, &w);
   vec_create(1, &wtmp);
   vec_create(v_size, &gtmp);
   vec_create(v_size, &out_v);
   vec_create(v_size, &out_w);
   vec_create(1, &out_vtmp);
   vec_create(1, &htmp);
   vec_create(v_size - 1, &out_u);
   double **app_u;
   for(int i = 0; i < ntime; i++)
   {
      
      //printf("\n\n\n %d\n", i);
      /* Compute Lu */
      app_u = app->u;
      wtmp[0] = 0.0;
      //printf("\n %f", wtmp[0]);
      //printf("\n %f", app_u[1000]);

      if (i != 0)
      {
         vec_axpy(1, -1, app_u[i-1], wtmp);
         //printf("\n what the heckin %f", app_u[i-1][0]);
      }
      if ( i != app->ntime - 1 )
      {
         vec_axpy(1, 1, app_u[i], wtmp);
         //printf("\n what the heckin %f", app_u[i][0]);
      }
      //printf("\n %f", wtmp[0]);

      /* wtmp = ( g - Lu )  */
      gtmp[0] = g[i];
      vec_axpy(1, -1, gtmp, wtmp);
      vec_scale(1, -1.0, wtmp);
      //printf("\n %f", wtmp[0]);
      /* wtmp = D^{-1} (g - Lu) */
      apply_Dinv(dt, wtmp, wtmp);
      vec_copy(1, wtmp, out_vtmp);
      //printf("\n %f", wtmp[0]);
      /* v is D^{-1} (g - Lu) */
      out_v[i] = -out_vtmp[0];

      /* D^{-T} V D^{-1} (g - Lu) */
      apply_V(dt, wtmp, u, v, i, app->ntime, u_in, u_out);
      apply_DAdjointinv(dt, wtmp, wtmp);

      /* D^{-T} h  */
      htmp[0] = h[i];
      apply_DAdjointinv(dt, htmp, htmp);
      vec_axpy(1, 1, htmp, wtmp);
      out_w[i] = -wtmp[0];
      if(i < v_size - 1)
      {
         out_u[i] = app_u[i][0];
      }
   }
   
   merge(e, out_u, out_v, out_w, v_size);

   for (int i = 0; i < (app->ntime); i++)
   {
      free(app->u[i]);
   }
   free(app->u);
   free(app);
   braid_Destroy(core);
}

/*Given  x, Computes Hx*/
void ApplyH(double *u, double *v, double dt, double *x, double *Hx, int n, double u_in, double u_out){
   double * utmp;
   double * vtmp;
   double * wtmp;
   double * rtmp;
   vec_create(1, &utmp);
   vec_create(1, &vtmp);
   vec_create(1, &wtmp);
   vec_create(1, &rtmp);

   double* uOfx = x;
   double* vOfx = x+(n-1);
   double* wOfx = x+(n+(n-1));

   int HxIndex=0;
   for ( int i = 0; i < n-1; i++ ) // U*u + B^T*v + L^T*w
   {
      rtmp[0] = 0.0;
      // U*u
      utmp[0] = 0.0;
      apply_U(dt, utmp, u, v, i);
      rtmp[0] = utmp[0]*uOfx[i];

      // + B^T*v
      vtmp[0] = 0.0;
      vtmp[0] -= dt*v[i]/(u[i]*u[i])*vOfx[i];
      vtmp[0] -= dt*v[i+1]/(u[i]*u[i])*vOfx[i+1];
      rtmp[0] += vtmp[0];

      // + L^T*w
      rtmp[0] += wOfx[i] - wOfx[i+1];
      Hx[HxIndex] = rtmp[0];
      ++HxIndex;
   }

   for ( int i = 0; i < n; i++ ) // B*u + V*v - D^T*w
   {
      rtmp[0] = 0.0;
      //B*u
      utmp[0] = 0.0;
      if ( i != n-1 ) utmp[0] -= dt*v[i]/(u[i]*u[i])*uOfx[i];
      if ( i != 0 ) utmp[0] -= dt*v[i]/(u[i-1]*u[i-1])*uOfx[i-1];
      rtmp[0] = utmp[0];

      // + V*v
      vtmp[0]=0.0;
      apply_V(dt,vtmp,u,v,i,n, u_in, u_out);
      rtmp[0] += vtmp[0]*vOfx[i];

      // - D^T*w
      wtmp[0]=wOfx[i];
      apply_DAdjoint(dt, wtmp, wtmp);
      vec_axpy(1, -1, wtmp, rtmp);
      Hx[HxIndex] = rtmp[0];
      ++HxIndex;
   }

   for ( int i = 0; i < n; i++ ) // L*u - D*v
   {
      rtmp[0] = 0.0;
      // L*u
      wtmp[0] = 0;
      if ( i != n-1 ) wtmp[0] += uOfx[i];
      if ( i != 0 ) wtmp[0] -= uOfx[i-1];
      rtmp[0] = wtmp[0];

      // - D*v
      vtmp[0] = vOfx[i];
      apply_D(dt,vtmp,vtmp);
      vec_axpy(1,-1,vtmp,rtmp);
      Hx[HxIndex] = rtmp[0];
      ++HxIndex;
   }

   vec_destroy(utmp);
   vec_destroy(vtmp);
   vec_destroy(wtmp);
   vec_destroy(rtmp);
}
   
double quick_norm(double *v, int n){
   double sum = 0.0;
   for (int i = 0; i < n; i++){
      sum += v[i]*v[i];
   }
   return sqrt(sum);
}
 

/* compute F(x)=\nabla F(u,v,w)*/
void applyF(double *u, double *v, double *w, double *g, double *Fx, double dt, int n){
   double *utmp, *vtmp, *wtmp; 
   
   /*Julia Code
       fu = zeros(nc-1)
    for i = 1:nc-1
        fu[i] = -(dt/2)*(v[i]^2 + v[i+1]^2)/(u[i]^2)  #(54)
    end
    fu = fu + L'*w
    */
   vec_create(n-1, &utmp);
   for (int i = 0; i < n-1; i++){
      utmp[i] = -(dt/2) * (v[i]*v[i] + v[i+1]*v[i+1]) / (u[i]*u[i]) + (w[i]-w[i+1]);
   }
   
   /* Julia Code
       # Apply KKT gradient-v block equations   #(55)
    fv = zeros(nc)
    fv[1] = dt*(1/uin + 1/u[1])*v[1]
    for i = 2:nc-1
        fv[i] = dt*(1/u[i-1] + 1/u[i])*v[i]
    end
    fv[nc] = dt*(1/u[nc-1] + 1/uout)*v[nc]
    fv = fv - D'*w
   */
   vec_create(n, &vtmp);
   for (int i = 0; i < n; i++){
      double ulm2 = ( i==0? 1.0/u_in : 1.0/u[i-1] );
      double urm2 = ( i==n-1? 1.0/u_out : 1.0/u[i] );
      vtmp[i] = dt * (ulm2 + urm2) * v[i] + dt*w[i];
   }

   /* Julia Code
    # Apply KKT gradient-w block equations
    fw = L*u - D*v
    fw[1]  = fw[1]  - uin
    fw[nc] = fw[nc] + uout
   */
   vec_create(n, &wtmp);
   for (int i = 0; i < n; i++){
      wtmp[i] = 0;
      if(i != n-1) wtmp[i] += u[i];
      if(i != 0) wtmp[i] -= u[i-1];
      wtmp[i] += dt*v[i]-g[i];
   }

   merge(Fx, utmp, vtmp, wtmp, n);
   vec_destroy(utmp);
   vec_destroy(vtmp);
   vec_destroy(wtmp);
}
 

/* Apply Jac[F(x)]^-1*/
void ApplyJacFxInv(double *Fx, double *u, double *v, double dt, int n, double u_in, double u_out){
   double *dx, *r, *e, *Hdx;
   vec_create(3*n-1, &dx);
   vec_create(3*n-1, &r);
   vec_create(3*n-1, &e);
   vec_create(3*n-1, &Hdx);

   for(int i = 0; i < 3 * n- 1; i++)
   {
      dx[i] = 0;
      e[i] = 0;
      Hdx[i] = 0;
      r[i] = 0;
   }

   double rtol = 1e-3;
   int max_itr = 2;
   //compute the first residual   
   vec_copy(3*n-1, Fx, r);
   ApplyH(u, v, dt, dx, Hdx, n, u_in, u_out);
   vec_axpy(3*n-1, -1.0, Hdx, r);

   double rnorm0 = quick_norm(r,3*n-1);
   double rnrel  = 1.0;
   int itr = 0;
   while (rnrel > rtol && itr < max_itr){
      /*one iteration carries out dx = dx + Hp^-1*(Fx-H*dx)*/
      // r <- Fx-H*dx
      vec_copy(3*n-1, Fx, r);
      ApplyH(u, v, dt, dx, Hdx, n, u_in, u_out);
      vec_axpy(3*n-1, -1.0, Hdx, r);  

      solveHp(e, r, u, v, dt, n, u_in, u_out);
      for ( int i = 0; i < 59; ++i )
      {
         printf("%d %d: % 1.14e\n", itr, i+1, e[i]);
      }
      // dx <- dx + e
      vec_axpy(3*n-1, 1.0, e, dx);
      itr++;
      rnrel = quick_norm(r,3*n-1) / rnorm0;
      printf("    linear iter = %d, rnrel = % 1.14e\n", itr, rnrel);
   }
   vec_copy(3*n-1, dx, Fx);
   vec_destroy(dx);
   vec_destroy(r);
   vec_destroy(e);
   vec_destroy(Hdx);
}

/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
   double dt;
   int ntime, max_newton_iter;
   int rank;
   int xsize;
   double * x, * y, *u, * v, *w, *g ;
  

   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

/* Define time domain */
   ntime  = 20;              /* Total number of time-steps */
   tstart = 0.0;             /* Beginning of time domain */
   tstop  = 1.0;             /* End of time domain*/
   xsize = 3 * ntime - 1;
/*  maximum iterations of Newton's method */
   max_newton_iter = 10;
   
   /* create vectors */
   vec_create(xsize, &x);
   vec_create(xsize, &y);
   vec_create(ntime - 1, &u);
   vec_create(ntime, &v);
   vec_create(ntime, &w);
   vec_create(ntime, &g);
   //initialize g
   for(int i = 0; i < ntime; i++)
   {
      g[i] = 0;
   }   
   //initialize y such that its norm is greater than residual tolerance
   for(int i = 0; i < xsize; i++)
   {
      y[i] = 1;
   }
   //initialize u in , u out, put them into g
   u_in = 1.0;
   u_out = 2.0;
   g[0] = u_in;
   g[ntime-1] -= u_out;

   //get dt
   dt = (tstop - tstart)/ntime;

   //initialize x vector to all 1's
   for(int i = 0; i < xsize; i++)
   {
      x[i]=1;
   }
   
   extract(x,  u, v, w, ntime);
   applyF(u, v, w, g, y, dt, ntime);
   double ynorm0 = quick_norm(y, xsize);
   double ynrel = 1.0;
   
   int iter = 0;
   while( ynrel >= 1.0e-6 && iter < max_newton_iter )
   {
      //apply Jacobian inverse
      //  x <- x - Jac(F(x))^{-1} * F(x)
      ApplyJacFxInv(y, u, v, dt, ntime, u_in, u_out);
      vec_axpy(xsize, -1.0, y, x);

      //apply F to x after extracting u, v, w, store in
      extract(x,  u, v, w, ntime);
      applyF(u, v, w, g, y, dt, ntime);
      ++iter;
      ynrel = quick_norm(y, xsize) / ynorm0;
      printf("iter = %d, fnrel = % 1.14e\n", iter, ynrel);
   }
   MPI_Finalize();
}

   //    braid_Core  core;
//  