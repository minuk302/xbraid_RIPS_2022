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
   double **u;           /* State vectors at each time point on my proc */
   double *diag_U;      /* current u  */
   double *diag_V;     /* current v */
   double *k_bar;         /* rhs */
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

void applyLAdjoint(int size, double *x, double *LTx){
    for (int i = 0; i < size-1; i++)
    {
        LTx[i] = x[i] - x[i+1];
    }
}

void vec_set(int size, double v, double *x){
    for (int i = 0; i < size; i++)
    {
        x[i]=v;
    }
}

void applyL(int size, double *x, double *Lx){
    vec_set(size, 0, Lx);
    for (int i = 0; i < size-1; i++)
    {
        Lx[i] += x[i];
        Lx[i+1] -= x[i];
    }
}

void applyDiag(int size, double *D, double *x){
    for (int i = 0; i < size; i++)
    {
        x[i] = D[i] * x[i];
    }
}

void applyB(int size, double *u, double *v, double *x, double *ret, double dt){
    vec_set(size, 0, ret);
    for (int i = 0; i < size-1; i++)
    {
        ret[i] += -dt*x[i]*v[i]/(u[i]*u[i]);
        ret[i+1] += -dt*x[i]*v[i+1]/(u[i]*u[i]);
    }
}

void applyBAdjoint(int size, double *u, double *v, double *x, double *ret, double dt){
    for (int i = 0; i < size-1; i++)
    {
        ret[i] = -dt * (v[i]*x[i] + v[i+1]*x[i+1]) / (u[i] * u[i]);
    }
}


/*--------------------------------------------------------------------------
 * TriMGRIT wrapper routines
 *--------------------------------------------------------------------------*/

/* Compute A(u) - \bar{k} - f */

int my_TriResidual(braid_App       app,
                   braid_Vector    uleft,
                   braid_Vector    uright,
                   braid_Vector    f,
                   braid_Vector    r,
                   braid_TriStatus status)
{
   double  t, tprev, tnext, dt;
   int     level, index;

   braid_TriStatusGetTriT(status, &t, &tprev, &tnext);
   braid_TriStatusGetLevel(status, &level);
   braid_TriStatusGetTIndex(status, &index);

   double *diag_U = (app->diag_U);
   double *diag_V = (app->diag_V);
   double *k_bar = (app->k_bar);
//    int n = (app->npoints);
   
   /* Get the time-step size */
   if (t < tnext)
   {
      dt = tnext - t;
   }
   else
   {
      dt = t - tprev;
   }

   double *rtmp;
   vec_create(1, &rtmp);

   /* Fake equation - no need to even compute anything here - just set residual to zero */
   if (index == 0)
   {
      rtmp[0] = 0.0;
      vec_copy(1, rtmp, (r->values));
      return 0;
   }
    
   /* Compute action of center block */
   double uc = (r->values)[0];
   double rez = (dt*diag_U[index] + 1 / dt *(diag_V[index] + diag_V[index+1]))*uc;

   /* Compute action of west block */
   if (index > 1)
   {
      double ul = (uleft->values)[0];
      rez -= 1 / dt * diag_V[index]*ul;
   }
   
   /* Compute action of east block */
   if (uright != NULL)
   {
      double ur = (uright->values)[0];
      rez -= 1 / dt * diag_V[index + 1]*ur;
   }
   
   /* subtract bar{k} from rtmp*/
   rez -= dt*k_bar[index];

   /*load rez into rtmp vector*/
   rtmp[0] = rez;
   
   /* Subtract rhs f */
   if (f!= NULL)
   {
      vec_axpy(1, -1.0, (f->values), rtmp);
   }

   /* Copy temporary residual vector into residual */
   vec_copy(1, rtmp, (r->values));
   /* Destroy temporary vectors */
   vec_destroy(rtmp);
   
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

   braid_TriStatusGetTriT(status, &t, &tprev, &tnext);

   /* Get the time-step size */
   if (t < tnext)
   {
      dt = tnext - t;
   }
   else
   {
      dt = t - tprev;
   }
   
   int index;
   braid_TriStatusGetTIndex(status, &index);
   
   /* Fake equation - no need to even compute anything here */
   if (index == 0)
   {
      return 0;
   }

   /* Create temporary vector */
   double *utmp;
   vec_create(1, &utmp);

   /* Initialize temporary solution vector */
   vec_copy(1, (u->values), utmp);
   
   /* Compute residual */
   my_TriResidual(app, uleft, uright, f, u, status);

   /* Apply center block preconditioner (multiply by C^-1) to -r
    *
    * Using C = v_{i}/dt + v_{i+1}/dt + dt*u_{i}
    *          
    */
   double * diag_U = app->diag_U;
   double * diag_V = app->diag_V;
   double scale = dt*diag_U[index] + 1 / dt * (diag_V[index] + diag_V[index+1]);

   /* Complete residual update */
   vec_axpy(1, -1.0/scale, (u->values), utmp);
   vec_copy(1, utmp, (u->values));

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
   u->values[0] = 1;
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

/*------------------------------------*/


/*Using Xbraid to solve the state-based Schur complement of Hp*/
//ntime = len(diag_U) = len(sol)
void solveHpHelper(int ntime, double *diag_U, double *diag_V, double *k_bar, int rank, double *sol)
{

    /*pad diag_U, diag_V, k_bar with 1 (ghost equation)*/
    double *pad_U, *pad_V, *pad_k;
    vec_create(ntime+1, &pad_U);
    vec_create(ntime+2, &pad_V);
    vec_create(ntime+1, &pad_k);
    pad_U[0] = pad_V[0] = pad_k[0] = 1;
    for (int i = 0; i < ntime; i++) pad_U[i+1] = diag_U[i];
    for (int i = 0; i < ntime+1; i++) pad_V[i+1] = diag_V[i];
    for (int i = 0; i < ntime+1; i++) pad_k[i+1] = k_bar[i];
    
    
   braid_Core  core;
   my_App     *app;

   /* Define time domain and step */
   double tstart = 0.0;             /* Beginning of time domain */
   double tstop  = 0.75;             /* End of time domain WARNING: NEED TO BE 1-dt in this problem!!*/
//   double dt     = (tstop-tstart)/ntime; 

   int         max_levels, min_coarse, nrelax, nrelaxc, cfactor, maxiter;
   int         access_level, print_level;
   double tol;
   
//    printf("begin H_p inner solve\n");

   /* Define some Braid parameters */
   max_levels     = 10;
   min_coarse     = 1;
   nrelax         = 1;
   nrelaxc        = 7;
   maxiter        = 50;
   cfactor        = 2;
   tol            = 1.0e-6;
   access_level   = 1;
   print_level    = 0;//was 2

   /* Set up the app structure */
   app = (my_App *) malloc(sizeof(my_App));
   app->myid     = rank;
   app->ntime    = ntime;
   app->u        = NULL;
   app->diag_U = pad_U;
   app->diag_V = pad_V;
   app->k_bar = pad_k;

   /* Initialize XBr-1aid */
   braid_InitTriMGRIT(MPI_COMM_WORLD, MPI_COMM_WORLD, tstart, tstop, ntime, app,
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

   /*print load the solution into sol*/
   for (int i = 1; i < ntime + 1; i++)
   {
      sol[i-1]=(app->u)[i][0];
   }

   free(app);
   braid_Destroy(core);

   vec_destroy(pad_U);
   vec_destroy(pad_V);
   vec_destroy(pad_k);
}

int rank;

void show_vec(int size, double *v){
    printf("[");
    for (int i = 0; i < size; i++)
    {  
       printf("%f ",v[i]); 
    }
    printf("]\n");
}

/*solves for e in H_p e = r*/
void solveHp(double *e, double *r, double *u, double *v, double dt, int v_size, double u_in, double u_out){
//    printf("begin solve H_p\n");

    vec_set(3*v_size-1, 0, e);
    /*compute the diagonal of D_{uu}J, D_{vv}J*/
   double *DuuJ, *DvvJ; 
   vec_create(v_size-1, &DuuJ);
   vec_create(v_size, &DvvJ);
   for (int i = 0; i < v_size - 1; i++){
        DuuJ[i] = dt*(v[i]*v[i] + v[i+1]*v[i+1]) / (u[i]*u[i]*u[i]);
   }
   for (int i = 0; i < v_size; i++){
        double ul = (i==0? u_in: u[i-1]);
        double ur = (i==v_size-1? u_out: u[i]);
        DvvJ[i] = dt*(1/ul + 1/ur);
   }
   
   
   /* extract k, h, g from r */
   double *k, *h, *g;
   vec_create(v_size-1, &k);
   vec_create(v_size, &h);
   vec_create(v_size, &g);
   extract(r, k, h, g, v_size);

   

   /*compute k bar, the rhs for the schur complement*/
    double *k_bar, *tmps, *tmpl;
   vec_create(v_size - 1, &k_bar);
   vec_create(v_size - 1, &tmps);
   vec_create(v_size, &tmpl);
   // k_bar <- k
   vec_copy(v_size - 1, k, k_bar);
   // k_bar += L^T D^{-T} h
    vec_copy(v_size, h, tmpl);
    vec_scale(v_size, -1/dt, tmpl);
    applyLAdjoint(v_size, tmpl, tmps);
    vec_axpy(v_size-1, 1, tmps, k_bar);
    // k_bar += L^T D^{-T} V D^{-1} g
    vec_copy(v_size, g, tmpl);
    vec_scale(v_size, 1/(dt*dt), tmpl);
    applyDiag(v_size, DvvJ, tmpl);
    applyLAdjoint(v_size, tmpl, tmps);
    vec_axpy(v_size-1, 1, tmps, k_bar);



    /* solve the schur complement*/
    double *sol_u;
    vec_create(v_size-1, &sol_u);
    solveHpHelper(v_size-1, DuuJ, DvvJ, k_bar, rank, sol_u);


    /* back subtitue sol_u to solve for sol_v and sol_w*/
    double *sol_v, *sol_w;
    vec_create(v_size, &sol_v);
    vec_create(v_size, &sol_w);
    //tmp <- L*u
    applyL(v_size, sol_u, tmpl);
    //sol_v <- g
    vec_copy(v_size, g, sol_v);
    //sol_v -= Lu
    vec_axpy(v_size, -1, tmpl, sol_v);
    //sol_v <- -D^{-1} sol_v
    vec_scale(v_size, 1/dt, sol_v);



   //sol_w <- D^{-T} V sol_v
    vec_copy(v_size, sol_v, sol_w);
    applyDiag(v_size,DvvJ,sol_w);
    vec_scale(v_size,-1/dt,sol_w);
    //sol_w <- sol_w - D^{-T} h
    vec_copy(v_size, h, tmpl);
    vec_scale(v_size, -1/dt, tmpl);
    vec_axpy(v_size, -1, tmpl, sol_w);
    
    merge(e, sol_u, sol_v, sol_w, v_size);

    vec_destroy(DuuJ);
    vec_destroy(DvvJ);
    vec_destroy(k);
    vec_destroy(h);
    vec_destroy(g);
    vec_destroy(k_bar);
    vec_destroy(tmps);
    vec_destroy(tmpl);
    vec_destroy(sol_u);
    vec_destroy(sol_v);
    vec_destroy(sol_w);
}

/*Given  x, Computes Hx*/
void ApplyH(double *u, double *v, double dt, double *x, double *Hx, int n, double u_in, double u_out){

    /*compute the diagonal of D_{uu}J, D_{vv}J*/
    double *DuuJ, *DvvJ; 
    vec_create(n-1, &DuuJ);
    vec_create(n, &DvvJ);
    for (int i = 0; i < n - 1; i++){
        DuuJ[i] = dt*(v[i]*v[i] + v[i+1]*v[i+1]) / (u[i]*u[i]*u[i]);
    }
    for (int i = 0; i < n; i++){
        double ul = (i==0? u_in: u[i-1]);
        double ur = (i==n-1? u_out: u[i]);
        DvvJ[i] = dt*(1/ul + 1/ur);
    }


    /*input vectors*/
    double *x_u, *x_v, *x_w;
    vec_create(n-1, &x_u);
    vec_create(n, &x_v);
    vec_create(n, &x_w);
    extract(x, x_u, x_v, x_w, n);
   
    /* tmp vectors*/
    double *utmp, *vtmp, *wtmp;
    vec_create(n, &utmp);//has one extra slot
    vec_create(n, &vtmp);
    vec_create(n, &wtmp);
   
    /*output vectors*/
    double *retu, *retv, *retw;
    vec_create(n-1, &retu);
    vec_create(n, &retv);
    vec_create(n, &retw);
   
   // D_{uu}J*u + B^T*v + L^T*w
    vec_copy(n-1, x_u, utmp);
    applyDiag(n-1, DuuJ, utmp);    
    applyBAdjoint(n, u, v, x_v, vtmp, dt);
    applyLAdjoint(n, x_w, wtmp);
    vec_copy(n-1, utmp, retu);
    vec_axpy(n-1, 1, vtmp, retu);
    vec_axpy(n-1, 1, wtmp, retu);

    // B*u + D_{vv}J*v - D^T*w
    applyB(n, u, v, x_u, utmp, dt);
    vec_copy(n, x_v, vtmp);
    applyDiag(n, DvvJ, vtmp);
    vec_copy(n, x_w, wtmp);
    vec_scale(n, dt, wtmp);
    vec_copy(n, utmp, retv);
    vec_axpy(n, 1, vtmp, retv);
    vec_axpy(n, 1, wtmp, retv);

    // Lu - Dv
    applyL(n, x_u, utmp);
    vec_copy(n, x_v, vtmp);
    vec_scale(n, dt, vtmp);
    vec_copy(n, utmp, retw);
    vec_axpy(n, 1, vtmp, retw);

    // Hx = [retu; retv; retw]
    merge(Hx, retu, retv, retw, n);

    vec_destroy(x_u);
    vec_destroy(x_v);
    vec_destroy(x_w);
    vec_destroy(utmp);
    vec_destroy(vtmp);
    vec_destroy(wtmp);
    vec_destroy(retu);
    vec_destroy(retv);
    vec_destroy(retw);
}
   
/**/
double quick_norm(double *v, int n){
   double sum = 0.0;
   for (int i = 0; i < n; i++){
      sum += v[i]*v[i];
   }
   return sqrt(sum);
}
 

/* compute F(x)=\nabla L(u,v,w)*/
void applyF(double *u, double *v, double *w, double *g, double *Fx, double dt, int n){
   double *utmp, *vtmp, *wtmp; 
   
   vec_create(n-1, &utmp);
   for (int i = 0; i < n-1; i++){
      utmp[i] = -(dt/2) * (v[i]*v[i] + v[i+1]*v[i+1]) / (u[i]*u[i]) + (w[i]-w[i+1]);
   }
   
   vec_create(n, &vtmp);
   for (int i = 0; i < n; i++){
      double ul = i==0? g[0] : u[i-1];
      double ur = i==n-1? -g[n-1] : u[i];
      vtmp[i] = dt * (1/ul + 1/ur) * v[i] + dt*w[i];
   }
   
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
 

/* Apply Jac[F(x)]^{-1} = H^{-1}*/
void ApplyJacFxInv(double *Fx, double *u, double *v, double dt, int n, double u_in, double u_out){
    printf("    Begin applying H^{-1}\n");
    printf("    rhs:\n    ");
    show_vec(11, Fx);
    printf("    u:\n    ");
    show_vec(3, u);
    printf("    v:\n    ");
    show_vec(4, v);

   double *dx, *r, *e, *Hdx;
   vec_create(3*n-1, &dx);
   vec_create(3*n-1, &r);
   vec_create(3*n-1, &e);
   vec_create(3*n-1, &Hdx);
    
    vec_set(3*n-1,0,dx);
   double rtol = 1e-6;
   int max_itr = 10;
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

      // e <- Hp^-1 * r     
      solveHp(e, r, u, v, dt, n, u_in, u_out);
      // dx <- dx + e
      vec_axpy(3*n-1, 1.0, e, dx);
      itr++;

      rnrel = quick_norm(r,3*n-1) / rnorm0;
      printf("        linear iter = %d, rnrel = % 1.14e\n", itr, rnrel);
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
   
   /* Initialize MPI */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
//    double * diag_U, * diag_V, *k_bar;
//    vec_create(4, &diag_U);
//    vec_create(5, &diag_V);
//    vec_create(4, &k_bar);
//    for (int i = 0; i < 4; i++){
//         diag_U[i] = 1;
//         diag_V[i] = 1;
//         k_bar[i] = 1;
//    }
//    diag_V[4] = 1;

   
//    solveHpHelper(diag_U, diag_V, k_bar, 5, rank);


    // int n = 4;
    // double u[3]={3, 5, 4};
    // double v[4]={1, 1, 1, 1};
    // double u_in = 1, u_out=2, dt=1.0/4;
    // double rhs[11]={1,2,3,4,5,6,7,8,9,10,11};
    
    // double *ans;
    // vec_create(11, &ans);
    // show_vec(11, rhs);
    // ApplyJacFxInv(rhs, u, v, dt, n, u_in, u_out);
    // show_vec(11, rhs);

    /* Define time domain */
    int ntime = 8;              /* Total number of time-steps */
    double tstart = 0.0;        /* Beginning of time domain */
    double tstop  = 1.0;        /* End of time domain*/
    double dt = (tstop - tstart)/ntime;
   
   
   int max_newton_iter = 20;   /*  maximum iterations of Newton's method */
   
   int xsize = 3 * ntime - 1;
   double * x, * xtmp, *u, * v, *w, *g;
  
   
   /* create vectors */
   vec_create(xsize, &x);
   vec_create(xsize, &xtmp);

   vec_create(ntime - 1, &u);
   vec_create(ntime, &v);
   vec_create(ntime, &w);

   vec_create(ntime, &g);
   

   //initialize u in , u out, put them into g
   double u_in = 1.0;
   double u_out = 2.0;
   vec_set(ntime, 0, g);
   g[0] = u_in;
   g[ntime-1] = -u_out;

   
   
   //a straight line feels like a good guess
   vec_set(xsize, 1, x);
   for(int i = 0; i < ntime; i++)
   {
      x[i] = 1 + i/ntime;
   }
   
   
   double ynorm0 = ntime;
   double ynrel = 1.0;
   
   int iter = 0;
   while(ynrel >= 1.0e-6 && iter < max_newton_iter)
   {
        // unpack x=[u; v; w]
        extract(x, u, v, w, ntime);
        // xtmp <- x
        vec_copy(xsize, x, xtmp);
        // xtmp <- F(x)
        applyF(u, v, w, g, xtmp, dt, ntime);
        printf("Fx:\n");
        show_vec(xsize, xtmp);
        // xtmp <- Jac(F(x))^{-1} * xtmp
        ApplyJacFxInv(xtmp, u, v, dt, ntime, u_in, u_out);
        //  x <- x - Jac(F(x))^{-1} * F(x)
        vec_axpy(xsize, -1.0, xtmp, x);
        //house keeping
        ++iter;
        ynrel = quick_norm(xtmp, xsize) / ynorm0;
        printf("Newton iter = %d, fnrel = % 1.14e\n", iter, ynrel);
   }

    extract(x, u, v, w, ntime);
    printf("u:\n");
    show_vec(ntime-1, u);
    printf("v:\n");
    show_vec(ntime, v);
    printf("w:\n");
    show_vec(ntime, w);
    MPI_Finalize();
}

   //    braid_Core  core;
//  