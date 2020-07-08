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

//
// Example:       ex-01-pp.cpp
//
// Interface:     C++
// 
// Requires:      C-language and C++ support     
//
// Compile with:  make ex-01-pp
//
// Help with:     ex-01-pp -help
//
// Sample run:    mpirun -np 2 ex-01-pp
//
// Description:   solve the scalar ODE 
//                   u' = lambda u, 
//                   with lambda=-1 and y(0) = 1
//
//                Same as ex-01, only implements more advanced XBraid features.
//                
//                When run with the default 10 time steps, the solution is:
//                $ ./ex-01-pp
//                $ cat ex-01.out.00*
//                  1.00000000000000e+00
//                  6.66666666666667e-01
//                  4.44444444444444e-01
//                  2.96296296296296e-01
//                  1.97530864197531e-01
//                  1.31687242798354e-01
//                  8.77914951989026e-02
//                  5.85276634659351e-02
//                  3.90184423106234e-02
//                  2.60122948737489e-02
//                  1.73415299158326e-02
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "braid.hpp"

// --------------------------------------------------------------------------
// User-defined routines and objects 
// --------------------------------------------------------------------------

// Define BraidVector, can contain anything, and be named anything
// --> Put all time-dependent information here
class MyBraidVector
{
public:
   // Each vector holds the scalar solution value at a particular time 
   double value;

   // Construct a BraidVector for a given double 
   MyBraidVector(double value_) : value(value_) { }

   // Deconstructor
   virtual ~MyBraidVector() {};

};

typedef MyBraidVector *BraidVector;


// Wrapper for BRAID's App object 
// --> Put all time INDEPENDENT information here
class MyBraidApp : public BraidApp<MyBraidVector>
{
protected:
   // BraidApp defines tstart, tstop, ntime and comm_t

public:
   // Constructor 
   MyBraidApp(MPI_Comm comm_t_, int rank_, double tstart_ = 0.0, double tstop_ = 1.0, int ntime_ = 100);
   
   // We will need the MPI Rank
   int rank;

   // Deconstructor
   virtual ~MyBraidApp() {};

   // Define all the Braid Wrapper routines
   // Note: braid_Vector == BraidVector*
   virtual int Step(BraidVector    u_,
                    BraidVector    ustop_,
                    BraidVector    fstop_,
                    BraidStepStatus &pstatus);

   virtual int Clone(BraidVector  u_,
                     BraidVector *v_ptr);

   virtual int Init(double        t,
                    BraidVector *u_ptr);

   virtual int Free(BraidVector u_);

   virtual int Sum(double       alpha,
                   BraidVector x_,
                   double       beta,
                   BraidVector y_);

   virtual int SpatialNorm(BraidVector  u_,
                           double       *norm_ptr);

   virtual int BufSize(int *size_ptr,
                       BraidBufferStatus  &status);

   virtual int BufPack(BraidVector  u_,
                       void         *buffer,
                       BraidBufferStatus  &status);

   virtual int BufUnpack(void         *buffer,
                         BraidVector *u_ptr,
                         BraidBufferStatus  &status);

   virtual int Access(BraidVector       u_,
                      BraidAccessStatus &astatus);

   // Not needed in this example
   virtual int Residual(BraidVector     u_,
                        BraidVector     r_,
                        BraidStepStatus &pstatus) { return 0; }

   // Not needed in this example
   virtual int Coarsen(BraidVector   fu_,
                       BraidVector  *cu_ptr,
                       BraidCoarsenRefStatus &status) { return 0; }

   // Not needed in this example
   virtual int Refine(BraidVector   cu_,
                      BraidVector  *fu_ptr,
                      BraidCoarsenRefStatus &status)  { return 0; }

};

// Braid App Constructor
MyBraidApp::MyBraidApp(MPI_Comm comm_t_, int rank_, double tstart_, double tstop_, int ntime_)
   : BraidApp(comm_t_, tstart_, tstop_, ntime_)
{
   rank = rank_;
}

// 
int MyBraidApp::Step(BraidVector    u,
                     BraidVector    ustop,
                     BraidVector    fstop,
                     BraidStepStatus &pstatus)
{
   double tstart;             // current time
   double tstop;              // evolve to this time

   // Get time step information
   pstatus.GetTstartTstop(&tstart, &tstop);

   // Use backward Euler to propagate solution
   (u->value) = 1./(1. + tstop-tstart)*(u->value);
   
   // no refinement
   pstatus.SetRFactor(1);
   
   return 0;

}

int MyBraidApp::Init(double        t,
                     BraidVector *u_ptr)
{
   BraidVector u = new MyBraidVector(0.0);
   if (t != tstart)
   {
      u->value = 0.456;
   }
   else
   {
      u->value = 1.0;
   }

   *u_ptr = u;
   return 0;

}

int MyBraidApp::Clone(BraidVector  u,
                      BraidVector *v_ptr)
{
   BraidVector v = new MyBraidVector(u->value); 
   *v_ptr = v;

   return 0;
}


int MyBraidApp::Free(BraidVector u)
{
   delete u;
   return 0;
}

int MyBraidApp::Sum(double       alpha,
                    BraidVector x,
                    double       beta,
                    BraidVector y)
{
   (y->value) = alpha*(x->value) + beta*(y->value);
   return 0;
}

int MyBraidApp::SpatialNorm(BraidVector  u,
                            double       *norm_ptr)
{
   double dot;
   dot = (u->value)*(u->value);
   *norm_ptr = sqrt(dot);
   return 0;
}

int MyBraidApp::BufSize(int                *size_ptr,
                        BraidBufferStatus  &status)                           
{
   *size_ptr = sizeof(double);
   return 0;
}

int MyBraidApp::BufPack(BraidVector       u,
                        void               *buffer,
                        BraidBufferStatus  &status)
{
   double *dbuffer = (double *) buffer;

   dbuffer[0] = (u->value);
   status.SetSize(sizeof(double));

   return 0;
}

int MyBraidApp::BufUnpack(void              *buffer,
                          BraidVector      *u_ptr,
                          BraidBufferStatus &status)
{
   double *dbuffer = (double *) buffer;
   
   BraidVector u = new MyBraidVector(dbuffer[0]); 
   *u_ptr = u;

   return 0;
}

int MyBraidApp::Access(BraidVector       u,
                       BraidAccessStatus &astatus)
{
   char       filename[255];
   FILE      *file;

   // Extract information from astatus
   int done, level, iter, index;
   double t;
   astatus.GetTILD(&t, &iter, &level, &done);
   astatus.GetTIndex(&index);
   
   // Print information to file
   sprintf(filename, "%s.%04d.%03d", "ex-01.out", index, rank);
   file = fopen(filename, "w");
   fprintf(file, "%.14e\n", (u->value));
   fflush(file);
   fclose(file);

   return 0;
}


// --------------------------------------------------------------------------
// Main driver
// --------------------------------------------------------------------------

int main (int argc, char *argv[])
{

   double        tstart, tstop;
   int           ntime, rank;

   // Define time domain: ntime intervals
   ntime  = 10;
   tstart = 0.0;
   tstop  = tstart + ntime/2.;
  
   // Initialize MPI
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   
   // set up app structure
   MyBraidApp app(MPI_COMM_WORLD, rank, tstart, tstop, ntime);

   // Initialize Braid Core Object and set some solver options
   BraidCore<MyBraidVector> core(MPI_COMM_WORLD, &app);
   core.SetPrintLevel(2);
   core.SetMaxLevels(2);
   core.SetAbsTol(1.0e-6);
   core.SetCFactor(-1, 2);
   
   // Run Simulation
   core.Drive();

   // Clean up
   MPI_Finalize();

   return (0);
}



