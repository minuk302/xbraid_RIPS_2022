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
   jm[0] = dt * (v[index] * v[index] + v[index+1] * v[index+1])/(u[index] *u[index] * u[index] );
}

/*------------------------------------*/
/* apply_V(dt, INSERT MATRIX HERE, u, v, INSERT INDEX HERE);*/
void
apply_V(double dt, double * rob, double *u, double *v, int index)
{
   if (index == 0)
   {
      rob[0] = dt * (1/ u_in + 1/u[0]);
   }
   if (index == v_size - 1)
   {
      rob[0] = dt * (1/ u_out + 1/u[index-1])
   }
   if (index != 0 && index != v_size - 1)
   {
      rob[0] = dt * (1/u[index-1]+1/u[index]);
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
   double *rtmp, *utmp, *htmp, *h_othertemp, *gtmp, *ktmp, *u, *v;
   int     level, index;

   u = app->input_u;
   v = app->input_v;
   
   braid_TriStatusGetTriT(status, &t, &tprev, &tnext);
   braid_TriStatusGetLevel(status, &level);
   braid_TriStatusGetTIndex(status, &index);

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

   /* rtmp = U_i u */
   vec_copy(1, (r->values), utmp);
   //apply_U(dt, utmp);
   apply_U(dt, utmp, u, v, index);
   vec_copy(1, utmp, rtmp);

   /* rtmp = rtmp + D_i^-T V_i D_i^-1 u */
   vec_copy(1, (r->values), utmp);
   apply_Dinv(dt, utmp, utmp);
   apply_V(dt, utmp, u, v, index);
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
      apply_V(dt, utmp, u, v, index);
      //apply_V(dt, gamma_1, gamma_2, utmp);
      apply_DAdjointinv(dt, utmp, utmp);
      apply_PhiAdjoint(dt, utmp);
      vec_axpy(1, 1.0, utmp, rtmp);
   }

   /* Compute action of west block */
   if (uleft != NULL)
   {
      /* rtmp = rtmp - D_i^{-T} V_i D_i^{-1} Phi_i uleft */
      vec_copy(1, (uleft->values), utmp);
      apply_Phi(dt, utmp);
      apply_Dinv(dt, utmp, utmp);
      apply_V(dt, utmp, u, v, index);
      //apply_V(dt, gamma_1, gamma_2, utmp);
      apply_DAdjointinv(dt, utmp, utmp);
      vec_axpy(1, -1.0, utmp, rtmp);
   }
   
   /* Compute action of east block */
   if (uright != NULL)
   {
      /* rtmp = rtmp - Phi_{i+1}^T D_{i+1}^{-T} V_{i+1} D_{i+1}^{-1} uright */
      vec_copy(1, (uright->values), utmp);
      apply_Dinv(dt, utmp, utmp);
      apply_V(dt, utmp, u, v, index);
      //apply_V(dt, gamma_1, gamma_2, utmp);
      apply_DAdjointinv(dt, utmp, utmp);
      apply_PhiAdjoint(dt, utmp);
      vec_axpy(1, -1.0, utmp, rtmp);
   }
   /* subtract rhs kbar (add -k - L^T D^{-T} V D^{-1} g - L^T D^{-T} h*/
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
   apply_V(dt, gtmp, gtmp);
   apply_DAdjointinv(dt, gtmp, gtmp);
   apply_Dinv(dt, g_othertemp, g_othertemp);
   apply_V(dt, g_othertemp, g_othertempp);
   apply_DAdjointinv(dt, g_othertemp, g_othertemp);
   vec_copy(1,g_othertemp,rtmp);
   vec_axpy(1, -1, gtmp, rtmp);

   vec_axpy(1, -1, ktemp, rtmp);
   
   
   if (f!= NULL)
   {
   vec_copy(1, rtmp, (r->values));
   }
   /* Destroy temporary vectors */
   vec_destroy(rtmp);
   vec_destroy(utmp);
   
   return 0;
}   

/*------------------------------------*/

void extract(double *r, double *k, double *h, double *g, int v_size) 
{
   for(int i = 0; i < v_size - 1; i++)
   {
      k[i] = r[i];
   }
   for(int i =0 i < v_size; i++)
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
   for(int i =0 i < v _size; i++)
   {
      r[i+v_size - 1] = h[i];
   }
   for(int i = 0; i < v_size; i++)
   {
      r[i + 2 * v_size - 1] = g[i];
   }
}

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
   double * input_u = app->input_u;
   double * input_v = app->input_v;
   if (uright != NULL)
   {
      rtmp[0] = input_v[index]/(dt * dt) + input_v[index+1]/(dt * dt) + input_u[index];
   }
   else
   {
      /* At the rightmost point, use a different center coefficient approximation */
      rtmp[0] = input_v[index]/(dt * dt) + input_u[index];
   }

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
   u->values[0] = 1;
   u->values[sizeof(my_Vector)] = 2;
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

//   {
//      char  filename[255];
//      FILE *file;
//      int  iter;
//      braid_AccessStatusGetIter(astatus, &iter);
//
//      braid_AccessStatusGetTIndex(astatus, &index);
//      sprintf(filename, "%s.%02d.%04d.%03d", "trischur-ex-04.out", iter, index, app->myid);
//      file = fopen(filename, "w");
//      fprintf(file, "%1.14e, %1.14e\n", (u->values)[0], (u->values)[1]);
//      fflush(file);
//      fclose(file);
//   }

   return 0;
}

/*------------------------------------*/

int
my_BufSize(braid_App           app,
           int                 *size_ptr,
           braid_BufferStatus  bstatus)
{
   *size_ptr = 2*sizeof(double);
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

   for(i = 0; i < 2; i++)
   {
      dbuffer[i] = (u->values)[i];
   }

   braid_BufferStatusSetSize( bstatus,  2*sizeof(double));

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
   u->values = (double*) malloc( 2*sizeof(double) );

   /* Unpack the buffer */
   for(i = 0; i < 2; i++)
   {
      (u->values)[i] = dbuffer[i];
   }

   *u_ptr = u;
   return 0;
}


/*solves for e in H_p e = r*/
void solveHp(double *e, double *r, double *u, double *v, double dt, int v_size){
   braid_Core  core;
   my_App     *app;
         
   double      tstart, tstop, dt; 
   int         rank, ntime, arg_index;
   int         max_levels, min_coarse, nrelax, nrelaxc, cfactor, maxiter;
   int         access_level, print_level;
   double      tol;

   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   /* Define time domain */
   ntime  = 20;              /* Total number of time-steps */
   tstart = 0.0;             /* Beginning of time domain */
   tstop  = 1.0;             /* End of time domain*/

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

   /* Solve adjoint equations starting at time point t1=dt, recognizing that
    * braid will label this time point as index=0 instead of 1 */
   dt = (tstop-tstart)/ntime;

   /* Set up twohe app structure */
   app = (my_App *) malloc(sizeof(my_App));
   app->myid     = rank;
   app->ntime    = ntime;
   app->u        = NULL;
   
   /* Initialize k, h, g vectors */
   double * k;
   double * h;
   double * g;
   /* Initialive XBr-1aid */
   vec_create(v_size-1, &k);
   vec_create(v_size, &h);
   vec_create(v_size, &g);
   extract(r, k, h, g, v_size);
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
   braid_SetCFactor(coreewton_, -1, cfactor);
   braid_SetAccessLevel(core, access_level);
   braid_SetPrintLevel( core, print_level);       
   braid_SetMaxIter(core, maxiter);
   braid_SetAbsTol(core, tol);

   /* Parallel-in-time TriMGRIT simulation */
   braid_Drive(core);

   

   /*   Compute control v and adjoint w from x */
   double * out_u;
   double * utmp;
   double * out_v, * out_vtmp, * w, * wtmp, *gtmp, *htmp;
   vec_create(app->npoints, &w);
   vec_create(1, &wtmp);
   vec_create(1, utmp);
   vec_create(v_size, &gtmp);
   vec_create(v_size, &out_v);
   vec_create(1, &out_vtmp);
   vec_create(1, &htmp);
   vec_create(v_size - 1, &out_u);

   for(int i = 0; i < app->npoints; i++)
   {
      /* Compute Lu */

      double **out_u = (app->u);
      if (i != 0)
      {
         vec_copy(1, out_u[i-1], utmp);
      }
      vec_copy(1, out_u[i], wtmp);
      vec_axpy(1, -1, utmp, wtmp);

      /* subtract from g bar  */
      vec_copy(1, g[i] , gtmp);
      vec_axpy(1, -1, wtmp, gtmp);

      /* D^{-1} (g - Lu) */
      apply_Dinv(dt, wtmp, wtmp);
      vec_copy(1, wtmp, out_vtmp);

      /* v is D^{-1} (g - Lu) */
      out_v[i] = -out_vtmp[0];

      /* D^{-T} V D^{-1} (g - Lu) */
      apply_V(dt, wtmp, u, v, i);
      apply_DAdjointinv(dt, wtmp, wtmp);

      /* D^{-T} h  */
      vec_copy(1, wtmp, out_vtmp);
      htmp[0] = h[i];
      apply_DAdjointinv(dt, htmp, htmp);
      vec_axpy(1, 1, htmp, wtmp);
      w[i] = -wtmp[0];
   }

   merge(e, out_u, out_v, w, v_size);
   






      /* Compute adjoint w from u and print to file */
      /* ZTODO: This requires communication to do correctly */
      // {
      //    double *u;
      //    sprintf(filename, "%s.%03d", "trischur-trasport-0d-state.out.u", (app->myid));
      //    file = fopen(filename, "w");
      //    vec_create(1, &u);
      //    for (i = 0; i < (app->npoints); i++)
      //    {
      //       double **u = (app->u);
      //       if ((i+1) < (app->npoints))
      //       {
      //          vec_copy(1, w[i+1], u);
      //          apply_PhiAdjoint(dt, u);
      //          vec_axpy(1, -1.0, w[i], u);
      //       }
      //       else
      //       {
      //          vec_copy(1, w[i], u);
      //          vec_scale(1, -1.0, u);
      //       }
      //       apply_Uinv(dt, u);

      //       index = (app->ilower) + i + 1;
      //       fprintf(file, "%05d: % 1.14e, % 1.14e\n", index, w[0]);
      //    }
      //    vec_destroy(u);
      //    fflush(file);
      //    fclose(file);
      // }

      /* Compute control v from state u and print to file */
      // {
      //    double *v;
      //    double * utmp;
      //    double **u = (app->u);
      //    vec_create(1,&utmp);

      //    sprintf(filename, "%s.%03d", "trischur-transport-0d-state.out.v", (app->myid));
      //    file = fopen(filename, "w");
      //    vec_create(1, &v);
      //    for (i = 0; i < (app->npoints); i++)
      //    {
      //       if (i > 0)
      //       {
      //          vec_copy(1, u[i-1], v);
      //          apply_Phi(dt, v);
      //          vec_scale(1, -1, v);
      //          vec_axpy(1, 1.0, u[i], v);
      //       }
      //       else
      //       {
      //          vec_copy(1, u[i], v);
      //       }
      //       //apply_Phi(dt, v);
      //       //vec_scale(1, -1.0, v);
      //         /* Subtract g */
      //       if (i == 0)
      //       {
      //          /* rtmp = rtmp + g; g = Phi_0 u_0 */   
      //          utmp[0] =  0.0;
      //          apply_Phi(dt, utmp);
      //          vec_axpy(1, -1.0, utmp, v);
      //       }
      //         if (i == (app-> npoints) - 1)
      //         {

      //         }
      //       apply_Dinv(dt, v, v);
      //       index = (app->ilower) + i + 1;
      //       fprintf(file, "%05d: % 1.14e, % 1.14e\n", index, v[0]);
      //    }
      //    vec_destroy(v);
         //    fflush(file);
      //    ust you :)lose(file);
//Minuk I trsu
      //    for (i = 0; i < (app->npoints); i++)
      //    {
      //       free(app->u[i]);
      //    }
      //    free(app->u);
      // }

   
}

/*    free(app);
    
   braid_Destroy(core);
   MPI_Finalize();

   return (0);
 
 }
*/

/*Given  x, Computes Hx*/
void ApplyH(double *u, double *v, double dt, double *x, double *Hx, int n){
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
      rtmp[0] = utmp[0]*uOfx[i]

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
      apply_V(dt,vtmp,u,v,i);
      rtmp[0] += vtmp[0]*vOfx[i]

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
      vec_axpy(1,-1,vtmp,rtmp)
      Hx[HxIndex] = rtmp[0];
      ++HxIndex;
   }

   vec_destroy(utmp);
   vec_destroy(vtmp);
   vec_destroy(wtmp);
   vec_destroy(rtmp);
}
   
void quick_norm(double *v, int n){
   double sum = 0.0;
   for (int i = 0; i < n; i++){
      sum += v[i]*v[i];
   }
   return sqrt(sum);
}
 

/* compute F(x)=\nabla F(u,v,w)*/
void ApplyF(double *u, double *v, double *w, double *g, double *Fx, double dt, int n){
   double *utmp, *vtmp, *wtmp; 
   
   vec_create(n-1, &utmp);
   for (int i = 0; i < n-1; i++){
      utmp[i] = -(dt/2) * (v[i]*v[i] + v[i+1]*v[i+1]) / (u[i]*u[i]) + (w[i]-w[i+1]);
   }
   vec_copy(n-1, utmp, u);
   
   vec_create(n, &vtmp);
   for (int i = 0; i < n; i++){
      double ulm2 = i==0? 1/(g[0]*g[0]) : 1/(u[i-1]*u[i-1]);
      double urm2 = i==n-1? 1/(g[n-1]*g[n-1]) : 1/(u[i]*u[i]);
      vtmp[i] = dt * (ulm2 + urm2) * v[i] + dt*w[i];
   }
   vec_copy(n, vtmp, v);

   vec_create(n, &wtmp);
   for (int i = 0; i < n; i++){
      wtmp[i] = 0;
      if(i != n-1) wtmp[i] += u[i];
      if(i != 0) wtmp[i] -= u[i-1];
      wtmp[i] += dt*v[i]-g[i];
   }
   vec_copy(n, wtmp, w);

   merge(Fx, utmp, vtmp, wtmp, n);
   vec_destroy(utmp);
   vec_destroy(vtmp);
   vec_destroy(wtmp);
}
 

/* Apply Jac[F(x)]^-1*/
void ApplyJacFxInv(double *Fx, double *u, double *v, double dt, int n){
   double *dx, *r, *e, *Hdx;
   vec_create(3*n-1, &dx);
   vec_create(3*n-1, &r);
   vec_create(3*n-1, &e);
   vec_create(3*n-1, &Hdx);

   double tol = 1e-6;
   int max_itr = 10;
   //compute the first residual   
   vec_copy(3*n-1, Fx, r);
   ApplyH(u, v, dt, dx, Hdx, n);
   vec_axpy(3*n-1, -1.0, Hdx, r);   
   int itr = 0;

   while (quick_norm(r,3*n-1) > tol && itr < max_itr){
      /*one iteration carries out dx = dx + Hp^-1*(Fx-H*dx)*/
      // r <- Fx-H*dx
      vec_copy(3*n-1, Fx, r);
      ApplyH(u, v, dt, dx, Hdx, n);
      vec_axpy(3*n-1, -1.0, Hdx, r);   
      // e <- Hp^-1 * r     
      solveHp(*e, *r, *u, *v, dt, n);
      // dx <- dx + e
      vec_axpy(3*n-1, 1.0, e, dx);
      itr++;
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
   double   tstart, tstop, dt;
   int ntime, max_newwton_iter;
   int xsize = 3 * ntime - 1;
   double * x, * y, *u, * v, *w, *g ;
   vec_create(xsize, x);
   vec_create(xsize, y);
   vec_create(ntime - 1, u);
   vec_create(ntime, v);
   vec_create(ntime, w);
   vec_create(ntime, g);

/* Define time domain */
    ntime  = 20;              /* Total number of time-steps */
    tstart = 0.0;             /* Beginning of time domain */
    tstop  = 1.0;             /* End of time domain*/

/*  maximum iterations of Newton's method */
   max_newton_iter = 10;
   
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
   g[ntime-1] = -u_out;

   //get dt
   dt = (tstop - tstart)/ntime;

   //initialize x vector to all 1's
   for(int i = 0; i < ntime; i++)
   {
      x[i]=1;
   }
   
   extract(x,  u, v, w, ntime);
   applyF(u, v, w, g, y, dt, ntime);
   double ynorm0 = quick_norm(xsize, y);
   double ynrel = 1.0;
   
   int iter = 0;
   while( nrel >= 1.0e-6 && iter < max_newwton_iter )
   {
      //apply Jacobian inverse
      //  x <- x - Jac(F(x))^{-1} * F(x)
      ApplyJacFxInv(y, u, v, dt, ntime);
      vec_axpy(xsize, -1.0, y, x);

      //apply F to x after extracting u, v, w, store in
      extract(x,  u, v, w, ntime);
      applyF(u, v, w, g, y, dt, ntime);

      ++iter;
      ynrel = quick_norm(xsize, y) / ynorm0;
      printf("iter = %d, fnrel = % 1.14e\n", iter, ynrel);
   }

   //    braid_Core  core;
//    my_App     *app;
//    double      tstart, tstop, dt; 
//    int         rank, ntime, arg_index;
//    int         max_levels, min_coarse, nrelax, nrelaxc, cfactor, maxiter, max_newton_iter;
//    int         access_level, print_level;
//    double      tol;

//    /* Initialize MPI */
//    MPI_Init(&argc, &argv);
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//    /* Define time domain */
//    ntime  = 20;              /* Total number of time-steps */
//    tstart = 0.0;             /* Beginning of time domain */
//    tstop  = 1.0;             /* End of time domain*/

//    /* Define some Braid parameters */
//    max_levels     = 30;
//    min_coarse     = 1;
//    nrelax         = 1;
//    nrelaxc        = 7;
//    max_newton_iter    = 20;
//    maxiter        = 20;
//    cfactor        = 2;
//    tol            = 1.0e-6;
//    access_level   = 1;
//    print_level    = 2;

//    /* Give some initial conditions */


//    /* Parse command line */
//    arg_index = 1;
//    while (arg_index < argc)
//    {
//       if ( strcmp(argv[arg_index], "-help") == 0 )
//       {
//          printf("\n");
//          printf(" Solves a simple optimal control problem in time-serial on [0, 1] \n\n");
//          printf("  min   \\int_0^1 u_1(t)^2 + u_2(t)^2 + gamma_1 c_1(t)^2 + gamma_2 c_2(t)^2 dt \n\n");
//          printf("  s.t.  d/dt u_1(t) = u_2(t) \n");
//          printf("        d/dt u_2(t) = -u_2(t) + c(t) \n\n");
//          printf("  -tstop <tstop>          : Upper integration limit for time\n");
//          printf("  -ntime <ntime>          : Num points in time\n");
//          printf("  -ml <max_levels>        : Max number of braid levels \n");
//          printf("  -nu  <nrelax>           : Num F-C relaxations\n");
//          printf("  -nuc <nrelaxc>          : Num F-C relaxations on coarsest grid\n");
//          printf("  -mi <max_newton_iter>           : Max iterations \n");
//          printf("  -cf <cfactor>           : Coarsening factor \n");
//          printf("  -tol <tol>              : Stopping tolerance \n");
//          printf("  -access <access_level>  : Braid access level \n");
//          printf("  -print <print_level>    : Braid print level \n");
//          printf("\n");
//          exit(1);
//       }
//       else if ( strcmp(argv[arg_index], "-ntime") == 0 )
//       {
//          arg_index++;
//          ntime = atoi(argv[arg_index++]);
//       }
//       else if ( strcmp(argv[arg_index], "-tstop") == 0 )
//       {
//          arg_index++;
//          tstop = atof(argv[arg_index++]);
//       }
//       else if ( strcmp(argv[arg_index], "-ml") == 0 )
//       {
//          arg_index++;
//          max_levels = atoi(argv[arg_index++]);
//       }
//       else if ( strcmp(argv[arg_index], "-nu") == 0 )
//       {
//          arg_index++;
//          nrelax = atoi(argv[arg_index++]);
//       }
//       else if ( strcmp(argv[arg_index], "-nuc") == 0 )
//       {
//          arg_index++;
//          nrelaxc = atoi(argv[arg_index++]);
//       }
//       else imp(argv[arg_index], "-mi") == 0 )
//       {
//        else if ( strcmp(argv[arg_index], "-mni") == 0 )
//       {
//          arg_index++;
//          max_newton_iter = atoi(argv[arg_index++]);
//       }
//          arg_index++;
//          max_newton_iter = atoi(argv[arg_index++]);
//       }
//       else if ( strcmp(argv[arg_index], "-cf") == 0 )
//       {
//          arg_index++;
//          cfactor = atoi(argv[arg_index++]);
//       }
//       else if ( strcmp(argv[arg_index], "-tol") == 0 )
//       {
//          arg_index++;
//          tol = atof(argv[arg_index++]);
//       }
//       else if ( strcmp(argv[arg_index], "-access") == 0 )
//       {
//          arg_index++;
//          access_level = atoi(argv[arg_index++]);
//       }
//       else if ( strcmp(argv[arg_index], "-print") == 0 )
//       {
//          arg_index++;
//          print_level = atoi(argv[arg_index++]);
//       }
//       else
//       {
//          printf("ABORTING: incorrect command line parameter %s\n", argv[arg_index]);
//          return (0);
//       }
//    }

//    /* Solve adjoint equations starting at time point t1=dt, recognizing that
//     * braid will label this time point as index=0 instead of 1 */
//    dt = (tstop-tstart)/ntime;

//    /* Set up twohe app structure */
//    app = (my_App *) malloc(sizeof(my_App));
//    app->myid     = rank;
//    app->ntime    = ntime;
//    app->u        = NULL;
   
//    /* Initialize XBraid */
//    braid_InitTriMGRIT(MPI_COMM_WORLD, MPI_COMM_WORLD, dt, tstop, ntime-1, app,
//                       my_TriResidual, my_TriSolve, my_Init, my_Clone, my_Free,
//                       my_Sum, my_SpatialNorm, my_Access,
//                       my_BufSize, my_BufPack, my_BufUnpack, &core);

//    /* Set some XBraid(_Adjoint) parameters */
//    braid_SetMaxLevels(core, max_levels);
//    braid_SetMinCoarse(core, min_coarse);
//    braid_SetNRelax(core, -1, nrelax);
//    if (max_levels > 1)
//    {
//       braid_SetNRelax(core, max_levels-1, nrelaxc); /* nrelax on coarsest level */
//    }
// Newton_n   braid_SetCFactor(coreewton_, -1, cfactor);
//    braid_SetAccessLevel(core, access_level);
//    braid_SetPrintLevel( core, print_level);       
//    braid_SetMaxIter(core, maxiter);
//    braid_SetAbsTol(core, tol);

//    /* Parallel-in-time TriMGRIT simulation */
//    braid_Drive(core);

//    if (access_level > 0)
//    {
//       char  filename[255];
//       FILE *file;
//       int   i, index;

//       /* Print state u to file */
//       {
//          sprintf(filename, "%s.%03d", "trischur-transport-0d-state.out.u", (app->myid));
//          file = fopen(filename, "w");
//          for (i = 0; i < (app->npoints); i++)
//          {
//             double **u = (app->u);

//             index = (app->ilower) + i + 1;
//             fprintf(file, "%05d: % 1.14e, % 1.14e\n", index, u[i][0], u[i][1]);
//          }
//          fflush(file);
//          fclose(file);
//       }

//       /* Compute adjoint w from u and print to file */
//       /* ZTODO: This requires communication to do correctly */
//       // {
//       //    double *u;

//       //    sprintf(filename, "%s.%03d", "trischur-trasport-0d-state.out.u", (app->myid));
//       //    file = fopen(filename, "w");
//       //    vec_create(1, &u);
//       //    for (i = 0; i < (app->npoints); i++)
//       //    {
//       //       double **u = (app->u);

//       //       if ((i+1) < (app->npoints))
//       //       {
//       //          vec_copy(1, w[i+1], u);
//       //          apply_PhiAdjoint(dt, u);
//       //          vec_axpy(1, -1.0, w[i], u);
//       //       }
//       //       else
//       //       {
//       //          vec_copy(1, w[i], u);
//       //          vec_scale(1, -1.0, u);
//       //       }
//       //       apply_Uinv(dt, u);

//       //       index = (app->ilower) + i + 1;
//       //       fprintf(file, "%05d: % 1.14e, % 1.14e\n", index, u[0], u[1]);
//       //    }
//       //    vec_destroy(u);
//       //    fflush(file);
//       //    fclose(file);
//       // }

//       /* Compute control v from state u and print to file */
//       // {
//       //    double *v;
//       //    double * utmp;
//       //    double **u = (app->u);
//       //    vec_create(1,&utmp);

//       //    sprintf(filename, "%s.%03d", "trischur-transport-0d-state.out.v", (app->myid));
//       //    file = fopen(filename, "w");
//       //    vec_create(1, &v);
//       //    for (i = 0; i < (app->npoints); i++)
//       //    {
//       //       if (i > 0)
//       //       {
//       //          vec_copy(1, u[i-1], v);
//       //          apply_Phi(dt, v);
//       //          vec_scale(1, -1, v);
//       //          vec_axpy(1, 1.0, u[i], v);
//       //       }
//       //       else
//       //       {
//       //          vec_copy(1, u[i], v);
//       //       }
//       //       //apply_Phi(dt, v);
//       //       //vec_scale(1, -1.0, v);
//       //         /* Subtract g */
//       //       if (i == 0)
//       //       {
//       //          /* rtmp = rtmp + g; g = Phi_0 u_0 */   
//       //          utmp[0] =  0.0;
//       //          utmp[1] = -1.0;
//       //          apply_Phi(dt, utmp);
//       //          vec_axpy(1, -1.0, utmp, v);
//       //       }
//       //       apply_Dinv(dt, v, v);
//       //       index = (app->ilower) + i + 1;
//       //       fprintf(file, "%05d: % 1.14e, % 1.14e\n", index, v[0], v[1]);
//       //    }
//       //    vec_destroy(v);
//       //    fflush(file);
//       //    fclose(file);

//       //    for (i = 0; i < (app->npoints); i++)
//       //    {
//       //       free(app->u[i]);
//       //    }
//       //    free(app->u);
//       // }

//    }

//    free(app);
   
//    braid_Destroy(core);
//    MPI_Finalize();

//    return (0);
}//