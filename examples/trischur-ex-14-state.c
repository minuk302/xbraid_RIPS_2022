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
   double   gamma_1;       /* Relaxation parameter for objective function */
   double   gamma_2;       /* Relaxation parameter for objective function */
   int      ntime;       /* Total number of time-steps (starting at time 0) */

   int      ilower;      /* Lower index for my proc */
   int      iupper;      /* Upper index for my proc */
   int      npoints;     /* Number of time points on my proc */
   double **u;           /* Adjoing vectors at each time point on my proc */

} my_App;


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
   u[0] = u[0] + dt*u[1];
   u[1] = u[1] - dt*u[1];
}

/*------------------------------------*/

void
apply_PhiAdjoint(double dt, double *w)
{
   w[1] = dt*w[0] + w[1] - dt*w[1];
}

/*------------------------------------*/

void
apply_U(double dt, double *u)
{
   u[0] *= 2*dt;
   u[1] *= 2*dt;
}

/*------------------------------------*/

void
apply_V(double dt, double gamma_1, double gamma_2, double *v)
{
   v[0] *= 2*gamma_1*dt;
   v[1] *= 2*gamma_2*dt;
}

/*------------------------------------*/

/* Maps v to w (which may be the same pointer) */
void
apply_Dinv(double dt, double *v, double *w)
{
   double vtmpone = v[0];
   double vtmptwo = v[1];
   w[0] = vtmpone/dt; 
   w[1] = vtmptwo/dt;
}

/*------------------------------------*/

/* Maps w to v (which may be the same pointer) */
void
apply_DAdjointinv(double dt, double *w, double *v)
{
   double wtmpone = w[0];
   double wtmptwo = w[1];
   v[0] = wtmpone/dt; 
   v[1] = wtmptwo/dt;
}
/* -------------------------------------------
 * TriMGRIT wrapper routines
 *--------------------------------------------------------------------------*/

/* Compute A(u) - f */

int
my_TriResidual(braid_App       app,
               braid_Vector    uleft,
               braid_Vector    uright,
               braid_Vector    f,
               braid_Vector    r,
               braid_TriStatus status)
{
   double  t, tprev, tnext, dt;
   double  gamma_1 = (app->gamma_1);
   double  gamma_2 = (app->gamma_2);
   double *rtmp, *utmp;
   int     level, index;
   
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
   vec_create(2, &rtmp);
   vec_create(2, &utmp);

   /* rtmp = U_i u */
   vec_copy(2, (r->values), utmp);
   apply_U(dt, utmp);
   vec_copy(2, utmp, rtmp);

   /* rtmp = rtmp + D_i^-T V_i D_i^-1 u */
   vec_copy(2, (r->values), utmp);
   apply_Dinv(dt, utmp, utmp);
   apply_V(dt, gamma_1, gamma_2, utmp);
   apply_DAdjointinv( dt, utmp, utmp);
   vec_axpy(2, 1.0, utmp, rtmp);

   /* rtmp = rtmp + Phi_{i+1}^T D_{i+1}^{-T} V_{i+1} D_{i+1}^{-1} Phi_{i+1} u */
   /* This term is zero at time 0, since Phi_0 = 0 */
   if (uright != NULL)
   {
      vec_copy(2, (r->values), utmp);
      apply_Phi(dt, utmp);
      apply_Dinv(dt, utmp, utmp);
      apply_V(dt, gamma_1, gamma_2, utmp);
      apply_DAdjointinv(dt, utmp, utmp);
      apply_PhiAdjoint(dt, utmp);
      vec_axpy(2, 1.0, utmp, rtmp);
   }

   /* Compute action of west block */
   if (uleft != NULL)
   {
      /* rtmp = rtmp - D_i^{-T} V_i D_i^{-1} Phi_i uleft */
      vec_copy(2, (uleft->values), utmp);
      apply_Phi(dt, utmp);
      apply_Dinv(dt, utmp, utmp);
      apply_V(dt, gamma_1, gamma_2, utmp);
      apply_DAdjointinv(dt, utmp, utmp);
      vec_axpy(2, -1.0, utmp, rtmp);
   }
   
   /* Compute action of east block */
   if (uright != NULL)
   {
      /* rtmp = rtmp - Phi_{i+1}^T D_{i+1}^{-T} V_{i+1} D_{i+1}^{-1} uright */
      vec_copy(2, (uright->values), utmp);
      apply_Dinv(dt, utmp, utmp);
      apply_V(dt, gamma_1, gamma_2, utmp);
      apply_DAdjointinv(dt, utmp, utmp);
      apply_PhiAdjoint(dt, utmp);
      vec_axpy(2, -1.0, utmp, rtmp);
   }
   /* subtract rhs kbar (add (L^T D^{-T} V D^{-1}) g)*/
   if (index == 0)
   {      
      utmp[0] = 0.0;
      utmp[1] = -1.0;
      apply_Phi(dt, utmp);
      apply_Dinv(dt, utmp, utmp);
      apply_V(dt, gamma_1, gamma_2, utmp);
      apply_DAdjointinv(dt, utmp, utmp);
      vec_axpy(2, -1.0, utmp, rtmp);
   }
   
   if (f!= NULL)
   {
      vec_axpy(2, -1.0, (f->values),rtmp);
   }
   
   /* Copy temporary residual vector into residual */
   vec_copy(2, rtmp, (r->values));
   
   /* Destroy temporary vectors */
   vec_destroy(rtmp);
   vec_destroy(utmp);
   
   return 0;
}   

/*------------------------------------*/

/* Solve A(u) = f */

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
   double gamma_1 = (app->gamma_1);
   double gamma_2 = (app->gamma_2);
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
   vec_create(2, &utmp);

   /* Initialize temporary solution vector */
   vec_copy(2, (u->values), utmp);
   
   /* Compute residual */
   my_TriResidual(app, uleft, uright, f, u, status);

   /* Apply center block preconditioner (multiply by \tilde{C}^-1) to -r
    *
    * Using \tilde{C} = | 4*gamma_1/dt + 2*dt            0             |
    *                   |  0                        4*gamma_2/dt + 2*dt |
    */
   rtmp = (u->values);
   if (uright != NULL)
   {
      rtmp[0] = -rtmp[0]/(4*gamma_1/dt + 2*dt);
      rtmp[1] = -rtmp[1]/(4*gamma_2/dt + 2*dt);
   }
   else
   {
      /* At the rightmost point, use a different center coefficient approximation */
      rtmp[0] = -rtmp[0]/(2*gamma_1/dt + 2*dt);
      rtmp[1] = -rtmp[1]/(2*gamma_2/dt + 2*dt);
   }

   /* Complete residual update */
   vec_axpy(2, 1.0, utmp, (u->values));
   
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
   vec_create(2, &(u->values));

   //u->values[0] = ((double)braid_Rand())/braid_RAND_MAX;
   //u->values[1] = ((double)braid_Rand())/braid_RAND_MAX;
   u->values[0] = 0.0;
   u->values[1] = -1.0;

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
   vec_create(2, &(v->values));

   /* Clone the values */
   v->values[0] = u->values[0];
   v->values[1] = u->values[1];

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
   (y->values)[1] = alpha*(x->values)[1] + beta*(y->values)[1];

   return 0;
}

/*------------------------------------*/

int
my_SpatialNorm(braid_App     app,
               braid_Vector  u,
               double       *norm_ptr)
{
   int i;
   double dot = 0.0;

   for (i = 0; i < 2; i++)
   {
      dot += (u->values)[i]*(u->values)[i];
   }
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
      vec_create(2, &(app->u[ii]));
      vec_copy(2, (u->values), (app->u[ii]));
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

/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/

int
main(int argc, char *argv[])
{
   braid_Core  core;
   my_App     *app;
         
   double      tstart, tstop, dt; 
   int         rank, ntime, arg_index;
   double      gamma_1, gamma_2;
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

   /* Define some optimization parameters */
   gamma_1 = 0.05;            /* Relaxation parameter in the objective function */
   gamma_2 = 0.05;            /* Relaxation parameter in the objective function */

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

   /* Parse command line */
   arg_index = 1;
   while (arg_index < argc)
   {
      if ( strcmp(argv[arg_index], "-help") == 0 )
      {
         printf("\n");
         printf(" Solves a simple optimal control problem in time-serial on [0, 1] \n\n");
         printf("  min   \\int_0^1 u_1(t)^2 + u_2(t)^2 + gamma_1 c_1(t)^2 + gamma_2 c_2(t)^2 dt \n\n");
         printf("  s.t.  d/dt u_1(t) = u_2(t) \n");
         printf("        d/dt u_2(t) = -u_2(t) + c(t) \n\n");
         printf("  -tstop <tstop>          : Upper integration limit for time\n");
         printf("  -ntime <ntime>          : Num points in time\n");
         printf("  -gamma_1 gamma_1>          : Relaxation parameter in the objective function \n");
         printf("  -gamma_2 <gamma_2>          : Relaxation parameter in the objective function \n");
         printf("  -ml <max_levels>        : Max number of braid levels \n");
         printf("  -nu  <nrelax>           : Num F-C relaxations\n");
         printf("  -nuc <nrelaxc>          : Num F-C relaxations on coarsest grid\n");
         printf("  -mi <maxiter>           : Max iterations \n");
         printf("  -cf <cfactor>           : Coarsening factor \n");
         printf("  -tol <tol>              : Stopping tolerance \n");
         printf("  -access <access_level>  : Braid access level \n");
         printf("  -print <print_level>    : Braid print level \n");
         printf("\n");
         exit(1);
      }
      else if ( strcmp(argv[arg_index], "-ntime") == 0 )
      {
         arg_index++;
         ntime = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-tstop") == 0 )
      {
         arg_index++;
         tstop = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-gamma_1") == 0 )
      {
         arg_index++;
         gamma_1 = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-gamma_2") == 0 )
      {
         arg_index++;
         gamma_2 = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-ml") == 0 )
      {
         arg_index++;
         max_levels = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nu") == 0 )
      {
         arg_index++;
         nrelax = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-nuc") == 0 )
      {
         arg_index++;
         nrelaxc = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-mi") == 0 )
      {
         arg_index++;
         maxiter = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-cf") == 0 )
      {
         arg_index++;
         cfactor = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-tol") == 0 )
      {
         arg_index++;
         tol = atof(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-access") == 0 )
      {
         arg_index++;
         access_level = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-print") == 0 )
      {
         arg_index++;
         print_level = atoi(argv[arg_index++]);
      }
      else
      {
         printf("ABORTING: incorrect command line parameter %s\n", argv[arg_index]);
         return (0);
      }
   }

   /* Solve adjoint equations starting at time point t1=dt, recognizing that
    * braid will label this time point as index=0 instead of 1 */
   dt = (tstop-tstart)/ntime;

   /* Set up twohe app structure */
   app = (my_App *) malloc(sizeof(my_App));
   app->myid     = rank;
   app->ntime    = ntime;
   app->gamma_1 = gamma_1;
   app->gamma_2 = gamma_2;
   app->u        = NULL;
   
   /* Initialize XBraid */
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

   if (access_level > 0)
   {
      char  filename[255];
      FILE *file;
      int   i, index;

      /* Print state u to file */
      {
         sprintf(filename, "%s.%03d", "trischur-ex-04-state.out.u", (app->myid));
         file = fopen(filename, "w");
         for (i = 0; i < (app->npoints); i++)
         {
            double **u = (app->u);

            index = (app->ilower) + i + 1;
            fprintf(file, "%05d: % 1.14e, % 1.14e\n", index, u[i][0], u[i][1]);
         }
         fflush(file);
         fclose(file);
      }

      /* Compute state u from adjoint w and print to file */
      /* ZTODO: This requires communication to do correctly */
      {
         double *w;
         double * utmp;
         double **u = (app->u);
         vec_create(2,&utmp);

         sprintf(filename, "%s.%03d", "trischur-ex-04-state.out.w", (app->myid));
         file = fopen(filename, "w");
         vec_create(2, &w);
         for (i = 0; i < (app->npoints); i++)
         {
            if (i > 0)
            {
               vec_copy(2, u[i-1], w);
               apply_Phi(dt, w);
               vec_scale(2, -1, w);
               vec_axpy(2, 1.0, u[i], w);
            }
            else
            {
               vec_copy(2, u[i], w);
            }
              /* Subtract g */
            if (i == 0)
            {
               /* rtmp = rtmp + g; g = Phi_0 u_0 */   
               utmp[0] =  0.0;
               utmp[1] = -1.0;
               apply_Phi(dt, utmp);
               vec_axpy(2, -1.0, utmp, w);
            }
            apply_Dinv(dt, w, w);
            apply_V(dt, gamma_1, gamma_2, w);
            apply_DAdjointinv(dt, w, w);
            index = (app->ilower) + i + 1;
            fprintf(file, "%05d: % 1.14e, % 1.14e\n", index, w[0], w[1]);
         }
         vec_destroy(w);
         fflush(file);
         fclose(file);
      }

      /* Compute control v from state u and print to file */
      {
         double *v;
         double * utmp;
         double **u = (app->u);
         vec_create(2,&utmp);

         sprintf(filename, "%s.%03d", "trischur-ex-04-state.out.v", (app->myid));
         file = fopen(filename, "w");
         vec_create(2, &v);
         for (i = 0; i < (app->npoints); i++)
         {
            if (i > 0)
            {
               vec_copy(2, u[i-1], v);
               apply_Phi(dt, v);
               vec_scale(2, -1, v);
               vec_axpy(2, 1.0, u[i], v);
            }
            else
            {
               vec_copy(2, u[i], v);
            }
              /* Subtract g */
            if (i == 0)
            {
               /* rtmp = rtmp + g; g = Phi_0 u_0 */   
               utmp[0] =  0.0;
               utmp[1] = -1.0;
               apply_Phi(dt, utmp);
               vec_axpy(2, -1.0, utmp, v);
            }
            apply_Dinv(dt, v, v);
            index = (app->ilower) + i + 1;
            fprintf(file, "%05d: % 1.14e, % 1.14e\n", index, v[0], v[1]);
         }
         vec_destroy(v);
         fflush(file);
         fclose(file);

         for (i = 0; i < (app->npoints); i++)
         {
            free(app->u[i]);
         }
         free(app->u);
      }

   }

   free(app);
   
   braid_Destroy(core);
   MPI_Finalize();

   return (0);
}