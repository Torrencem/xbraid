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
 * Example:       trischur-ex-04.cpp
 *
 * Interface:     C++
 * 
 * Requires:      C-language and C++ support
 *
 * Compile with:  make trischur-ex-04-pp
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

#include "tribraid.hpp"
#include "braid_test.h"

/* Define the state vector at one time-step */
class BraidVector
{
public:
    double *values; /* Holds the R^2 state vector (u_1, u_2) */

    BraidVector(double *values_) : values(values_) { }

    virtual ~BraidVector() {};
};

class MyBraidApp: public TriBraidApp
{
protected:
    // TriBraidApp defines tstart, tstop

public:
    // Constructor
    MyBraidApp(MPI_Comm comm_t_, int rank_, double gamma_, double tstart_ = 0.0, double tstop_ = 1.0, int ntime_ = 100);

    int myid;             /* Rank of the processor */
    double   gamma;       /* Relaxation parameter for objective function */
    int      ntime;       /* Total number of time-steps (starting at time 0) */

    int      ilower;      /* Lower index for my proc */
    int      iupper;      /* Upper index for my proc */
    int      npoints;     /* Number of time points on my proc */
    double **w;           /* Adjoint vectors at each time point on my proc */

    // Deconstructor
    virtual ~MyBraidApp() {};

    // Define all the TriBraid Wrapper routines
    // Note: braid_Vector == BraidVector*
    virtual int TriResidual(braid_Vector    uleft,
                            braid_Vector    uright,
                            braid_Vector    f,
                            braid_Vector    r,
                            BraidTriStatus  &status);

    virtual int TriSolve(braid_Vector    uleft,
                         braid_Vector    uright,
                         braid_Vector    fleft,
                         braid_Vector    fright,
                         braid_Vector    f,
                         braid_Vector    u,
                         BraidTriStatus  &status);

   virtual int Clone(braid_Vector  u_,
                     braid_Vector *v_ptr);

   virtual int Init(double        t,
                    braid_Vector *u_ptr);

   virtual int Free(braid_Vector u_);

   virtual int Sum(double       alpha,
                   braid_Vector x_,
                   double       beta,
                   braid_Vector y_);

   virtual int SpatialNorm(braid_Vector  u_,
                           double       *norm_ptr);

   virtual int BufSize(int *size_ptr,
                       BraidBufferStatus  &status);

   virtual int BufPack(braid_Vector  u_,
                       void         *buffer,
                       BraidBufferStatus  &status);

   virtual int BufUnpack(void         *buffer,
                         braid_Vector *u_ptr,
                         BraidBufferStatus  &status);

   virtual int Access(braid_Vector       u_,
                      BraidAccessStatus &astatus);

   // Not needed
   virtual int Residual(braid_Vector     u_,
                        braid_Vector     r_,
                        BraidStepStatus &pstatus) { return 0; }
};

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
   w[1] = w[1] - dt*w[1] + dt*w[0];
}

/*------------------------------------*/

void
apply_Uinv(double dt, double *u)
{
   u[0] /= 2*dt;
   u[1] /= 2*dt;
}

/*------------------------------------*/

void
apply_Vinv(double dt, double gamma, double *v)
{
   v[0] /= 2*gamma*dt;
}

/*------------------------------------*/

/* Maps v to w (which may be the same pointer) */
void
apply_D(double dt, double *v, double *w)
{
   double vtmp = v[0];
   w[0] = 0.0;
   w[1] = -dt*vtmp;
}

/*------------------------------------*/

/* Maps w to v (which may be the same pointer) */
void
apply_DAdjoint(double dt, double *w, double *v)
{
   v[0] = -dt*w[1];
}

/*--------------------------------------------------------------------------
 * TriMGRIT wrapper routines
 *--------------------------------------------------------------------------*/

MyBraidApp::MyBraidApp(MPI_Comm comm_t_, int rank_, double gamma_, double tstart_, double tstop_, int ntime_) : TriBraidApp(comm_t_, tstart_, tstop_, ntime_) {
    myid = rank_;
    gamma = gamma_;
    tstart = tstart_;
    tstop = tstop_;
    ntime = ntime_;
}

/* Compute A(u) - f */

int
MyBraidApp::TriResidual(braid_Vector    uleft_,
                        braid_Vector    uright_,
                        braid_Vector    f_,
                        braid_Vector    r_,
                        BraidTriStatus  &status)
{
   BraidVector *uleft = (BraidVector*) uleft_;
   BraidVector *uright = (BraidVector*) uright_;
   BraidVector *f = (BraidVector*) f_;
   BraidVector *r = (BraidVector*) r_;
   double  t, tprev, tnext, dt;
   double *rtmp, *utmp;
   int     level, index;
   
   status.GetTriT(&t, &tprev, &tnext);
   status.GetLevel(&level);
   status.GetTIndex(&index);

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

   /* Compute action of center block */

   /* rtmp = U_i^{-1} u */
   vec_copy(2, (r->values), utmp);
   apply_Uinv(dt, utmp);
   vec_copy(2, utmp, rtmp);

   /* rtmp = rtmp + D_i^T V_i^{-1} D_i^T u */
   vec_copy(2, (r->values), utmp);
   apply_DAdjoint(dt, utmp, utmp);
   apply_Vinv(dt, gamma, utmp);
   apply_D(dt, utmp, utmp);
   vec_axpy(2, 1.0, utmp, rtmp);

   /* rtmp = rtmp + Phi_i U_{i-1}^{-1} Phi_i^T u */
   /* This term is zero at time 0, since Phi_0 = 0 */
   if (uleft != NULL)
   {
      vec_copy(2, (r->values), utmp);
      apply_PhiAdjoint(dt, utmp);
      apply_Uinv(dt, utmp);
      apply_Phi(dt, utmp);
      vec_axpy(2, 1.0, utmp, rtmp);
   }

   /* Compute action of west block */
   if (uleft != NULL)
   {
      /* rtmp = rtmp - Phi_i U_{i-1}^{-1} uleft */
      vec_copy(2, (uleft->values), utmp);
      apply_Uinv(dt, utmp);
      apply_Phi(dt, utmp);
      vec_axpy(2, -1.0, utmp, rtmp);
   }
   
   /* Compute action of east block */
   if (uright != NULL)
   {
      /* rtmp = rtmp - U_i^{-1} Phi_{i+1}^T uright */
      vec_copy(2, (uright->values), utmp);
      apply_PhiAdjoint(dt, utmp);
      apply_Uinv(dt, utmp);
      vec_axpy(2, -1.0, utmp, rtmp);
   }

   /* Subtract rhs gbar (add g) */
   if (index == 0)
   {
      /* rtmp = rtmp + g; g = Phi_0 u_0 */
      utmp[0] =  0.0;
      utmp[1] = -1.0;
      apply_Phi(dt, utmp);
      vec_axpy(2, 1.0, utmp, rtmp);
   }

   /* Subtract rhs f */
   if (f != NULL)
   {
      /* rtmp = rtmp - f */
      vec_axpy(2, -1.0, (f->values), rtmp);
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
MyBraidApp::TriSolve(braid_Vector    uleft_,
                     braid_Vector    uright_,
                     braid_Vector    fleft_,
                     braid_Vector    fright_,
                     braid_Vector    f_,
                     braid_Vector    u_,
                     BraidTriStatus &status)
{
   BraidVector *uleft = (BraidVector *)uleft_;
   BraidVector *uright = (BraidVector *)uright_;
   // BraidVector *fleft = (BraidVector *)fleft_;
   // BraidVector *fright = (BraidVector *)fright_;
   BraidVector *f = (BraidVector *)f_;
   BraidVector *u = (BraidVector *)u_;
   double  t, tprev, tnext, dt;
   double *utmp, *rtmp;
   
   /* Get the time-step size */
   status.GetTriT(&t, &tprev, &tnext);
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
   TriResidual((braid_Vector) uleft,(braid_Vector) uright, (braid_Vector) f, (braid_Vector) u, status);

   /* Apply center block preconditioner (multiply by \tilde{C}^-1) to -r
    *
    * Using \tilde{C} = | 1/dt            0             |
    *                   |  0    ( 1/dt + dt/(2*gamma) ) |
    */
   rtmp = (u->values);
   if (uleft != NULL)
   {
      rtmp[0] = -rtmp[0]*dt;
      rtmp[1] = -rtmp[1]/(1/dt + dt/(2*gamma));
   }
   else
   {
      /* At the leftmost point, use a different center coefficient approximation */
      rtmp[0] = -rtmp[0]*(2*dt);
      rtmp[1] = -rtmp[1]/(1/(2*dt) + dt/(2*gamma));
   }

   /* Complete residual update */
   vec_axpy(2, 1.0, utmp, (u->values));
   
   /* no refinement */
   status.SetRFactor(1);

   /* Destroy temporary vectors */
   vec_destroy(utmp);

   return 0;
}   

/*------------------------------------*/

/* This is only called from level 0 */

int
MyBraidApp::Init(double        t,
                 braid_Vector *u_ptr)
{
   BraidVector *u;

   /* Allocate the vector */
   u = (BraidVector *) malloc(sizeof(BraidVector));
   vec_create(2, &(u->values));

   u->values[0] = ((double)braid_Rand())/braid_RAND_MAX;
   u->values[1] = ((double)braid_Rand())/braid_RAND_MAX;
//   u->values[0] = 1.0;
//   u->values[1] = 0.0;

   *u_ptr = (braid_Vector) u;

   return 0;
}

/*------------------------------------*/

int
MyBraidApp::Clone(braid_Vector  u_,
                  braid_Vector *v_ptr)
{
   BraidVector *u = (BraidVector *) u_;
   BraidVector *v;

   /* Allocate the vector */
   v = (BraidVector *) malloc(sizeof(BraidVector));
   vec_create(2, &(v->values));

   /* Clone the values */
   v->values[0] = u->values[0];
   v->values[1] = u->values[1];

   *v_ptr = (braid_Vector) v;

   return 0;
}

/*------------------------------------*/

int
MyBraidApp::Free(braid_Vector u_)
{
   BraidVector *u = (BraidVector *) u_;
   free(u->values);
   free(u);

   return 0;
}

/*------------------------------------*/

int
MyBraidApp::Sum(double        alpha,
                braid_Vector  x_,
                double        beta,
                braid_Vector  y_)
{
   BraidVector *x = (BraidVector *) x_;
   BraidVector *y = (BraidVector *) y_;
   (y->values)[0] = alpha*(x->values)[0] + beta*(y->values)[0];
   (y->values)[1] = alpha*(x->values)[1] + beta*(y->values)[1];

   return 0;
}

/*------------------------------------*/

int
MyBraidApp::SpatialNorm(braid_Vector  u_,
                        double       *norm_ptr)
{
   BraidVector *u = (BraidVector *) u_;
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
MyBraidApp::Access(braid_Vector       u_,
                   BraidAccessStatus &astatus)
{
   BraidVector *u = (BraidVector *) u_;
   int   done, index, ii;

   /* Print solution to file if simulation is over */
   astatus.GetDone(&done);

   if (done)
   {
      astatus.GetILowerUpper(&(this->ilower), &(this->iupper));
      (this->npoints) = (this->iupper) - (this->ilower) + 1;

      /* Allocate w array in app */
      if ((this->w) == NULL)
      {
         (this->w) = (double **) calloc((this->npoints), sizeof(double *));
      }

      astatus.GetTIndex(&index);
      ii = index - (this->ilower);
      if (this->w[ii] != NULL)
      {
         free(this->w[ii]);
      }
      vec_create(2, &(this->w[ii]));
      vec_copy(2, (u->values), (this->w[ii]));
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
MyBraidApp::BufSize(int                 *size_ptr,
                    BraidBufferStatus  &bstatus)
{
   *size_ptr = 2*sizeof(double);
   return 0;
}

/*------------------------------------*/

int
MyBraidApp::BufPack(braid_Vector        u_,
                    void               *buffer,
                    BraidBufferStatus  &bstatus)
{
   BraidVector *u = (BraidVector *) u_;
   double *dbuffer = (double *) buffer;
   int i;

   for(i = 0; i < 2; i++)
   {
      dbuffer[i] = (u->values)[i];
   }

   bstatus.SetSize( 2*sizeof(double));

   return 0;
}

/*------------------------------------*/

int
MyBraidApp::BufUnpack(void               *buffer,
                      braid_Vector       *u_ptr,
                      BraidBufferStatus  &bstatus)
{
   BraidVector *u = NULL;
   double    *dbuffer = (double *) buffer;
   int i;

   /* Allocate memory */
   u = (BraidVector *) malloc(sizeof(BraidVector));
   u->values = (double*) malloc( 2*sizeof(double) );

   /* Unpack the buffer */
   for(i = 0; i < 2; i++)
   {
      (u->values)[i] = dbuffer[i];
   }

   *u_ptr = (braid_Vector) u;
   return 0;
}

/*--------------------------------------------------------------------------
 * Main driver
 *--------------------------------------------------------------------------*/

int
main(int argc, char *argv[])
{
   double      tstart, tstop, dt; 
   int         rank, ntime, arg_index;
   double      gamma;
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
   gamma = 0.005;            /* Relaxation parameter in the objective function */

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
         printf("  min   \\int_0^1 u_1(t)^2 + u_2(t)^2 + gamma c(t)^2  dt \n\n");
         printf("  s.t.  d/dt u_1(t) = u_2(t) \n");
         printf("        d/dt u_2(t) = -u_2(t) + c(t) \n\n");
         printf("  -tstop <tstop>          : Upper integration limit for time\n");
         printf("  -ntime <ntime>          : Num points in time\n");
         printf("  -gamma <gamma>          : Relaxation parameter in the objective function \n");
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
         tstop = atoi(argv[arg_index++]);
      }
      else if ( strcmp(argv[arg_index], "-gamma") == 0 )
      {
         arg_index++;
         gamma = atof(argv[arg_index++]);
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

   /* Set up the app structure */
   MyBraidApp app (MPI_COMM_WORLD, rank, gamma, tstart, tstop, ntime);
   app.myid     = rank;
   app.ntime    = ntime;
   app.gamma    = gamma;
   app.w        = NULL;
    
   BraidTriCore core(MPI_COMM_WORLD, &app);

   /* Set some XBraid(_Adjoint) parameters */
   core.SetMaxLevels(max_levels);
   core.SetMinCoarse(min_coarse);
   core.SetNRelax(-1, nrelax);
   if (max_levels > 1)
   {
      core.SetNRelax(max_levels-1, nrelaxc); /* nrelax on coarsest level */
   }
   core.SetCFactor(-1, cfactor);
   core.SetAccessLevel(access_level);
   core.SetPrintLevel(print_level);       
   core.SetMaxIter(maxiter);
   core.SetAbsTol(tol);

   /* Parallel-in-time TriMGRIT simulation */
   core.Drive();

   if (access_level > 0)
   {
      char  filename[255];
      FILE *file;
      int   i, index;

      /* Print adjoint w to file */
      {
         sprintf(filename, "%s.%03d", "trischur-ex-04.out.w", (app.myid));
         file = fopen(filename, "w");
         for (i = 0; i < (app.npoints); i++)
         {
            double **w = (app.w);

            index = (app.ilower) + i + 1;
            fprintf(file, "%05d: % 1.14e, % 1.14e\n", index, w[i][0], w[i][1]);
         }
         fflush(file);
         fclose(file);
      }

      /* Compute state u from adjoint w and print to file */
      /* ZTODO: This requires communication to do correctly */
      {
         double *u;

         sprintf(filename, "%s.%03d", "trischur-ex-04.out.u", (app.myid));
         file = fopen(filename, "w");
         vec_create(2, &u);
         for (i = 0; i < (app.npoints); i++)
         {
            double **w = (app.w);

            if ((i+1) < (app.npoints))
            {
               vec_copy(2, w[i+1], u);
               apply_PhiAdjoint(dt, u);
               vec_axpy(2, -1.0, w[i], u);
            }
            else
            {
               vec_copy(2, w[i], u);
               vec_scale(2, -1.0, u);
            }
            apply_Uinv(dt, u);

            index = (app.ilower) + i + 1;
            fprintf(file, "%05d: % 1.14e, % 1.14e\n", index, u[0], u[1]);
         }
         vec_destroy(u);
         fflush(file);
         fclose(file);
      }

      /* Compute control v from adjoint w and print to file */
      {
         double *v;

         sprintf(filename, "%s.%03d", "trischur-ex-04.out.v", (app.myid));
         file = fopen(filename, "w");
         vec_create(2, &v);
         for (i = 0; i < (app.npoints); i++)
         {
            double **w = (app.w);

            apply_DAdjoint(dt, w[i], v);
            vec_scale(1, -1.0, v);
            apply_Vinv(dt, (app.gamma), v);

            index = (app.ilower) + i + 1;
            fprintf(file, "%05d: % 1.14e\n", index, v[0]);
         }
         vec_destroy(v);
         fflush(file);
         fclose(file);

         for (i = 0; i < (app.npoints); i++)
         {
            free(app.w[i]);
         }
         free(app.w);
      }

   }
   
   MPI_Finalize();

   return (0);
}
