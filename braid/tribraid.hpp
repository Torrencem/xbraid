// Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// Produced at the Lawrence Livermore National Laboratory.
// 
// This file is part of XBraid. For support, post issues to the XBraid Github page.
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License (as published by the Free Software
// Foundation) version 2.1 dated February 1999.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE. See the terms and conditions of the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public License along
// with this program; if not, write to the Free Software Foundation, Inc., 59
// Temple Place, Suite 330, Boston, MA 02111-1307 USA
//

#ifndef tribraid_hpp_HEADER
#define tribraid_hpp_HEADER

#include "_braid.h"
#include "braid.h"
#include "braid_status.h"
#include "braid_test.h"
#include "braid.hpp"

class BraidTriStatus;

// Wrapper for BRAID's App object. Users should inherit this class and implement
// the purely virtual functions (see braid.h for descriptions).
class TriBraidApp {
public:
   // User must set these four values as they are used by Braid.
   // This is done by the BraidApp constructor.
   MPI_Comm     comm_t;
   braid_Real   tstart;
   braid_Real   tstop;
   braid_Int    ntime;

   TriBraidApp(MPI_Comm _comm_t, braid_Real _tstart = 0.0,
               braid_Real _tstop = 1.0, braid_Int _ntime = 100)
       : comm_t(_comm_t), tstart(_tstart), tstop(_tstop), ntime(_ntime) {}

   virtual ~TriBraidApp() {}

   virtual braid_Int TriResidual(braid_Vector uleft_,
                                 braid_Vector uright_,
                                 braid_Vector f_,
                                 braid_Vector r_,
                                 BraidTriStatus &pstatus) = 0;

   virtual braid_Int TriSolve(braid_Vector uleft_,
                              braid_Vector uright_,
                              braid_Vector fleft_,
                              braid_Vector fright_,
                              braid_Vector f_,
                              braid_Vector u_,
                              BraidTriStatus &pstatus) = 0;

   /** @brief Compute the residual at time @a tstop, given the approximate
       solutions at @a tstart and @a tstop. The values of @a tstart and @a tstop
       can be obtained from @a pstatus.

       @param[in]     u_ Input: approximate solution at time @a tstop.
       @param[in,out] r_ Input: approximate solution at time @a tstart.
                         Output: residual at time @a tstop.

       @see braid_PtFcnResidual.
   */
   virtual braid_Int Residual(braid_Vector     u_,
                              braid_Vector     r_,
                              BraidStepStatus &pstatus) = 0;

   /// Allocate a new vector in @a *v_ptr, which is a deep copy of @a u_.
   /// @see braid_PtFcnClone.
   virtual braid_Int Clone(braid_Vector  u_,
                           braid_Vector *v_ptr) = 0;

   /** @brief Allocate a new vector in @a *u_ptr and initialize it with an
       initial guess appropriate for time @a t. If @a t is the starting time,
       this method should set @a *u_ptr to the initial value vector of the ODE
       problem.
       @see braid_PtFcnInit. */
   virtual braid_Int Init(braid_Real    t,
                          braid_Vector *u_ptr) = 0;

   /// De-allocate the vector @a u_.
   /// @see braid_PtFcnFree.
   virtual braid_Int Free(braid_Vector u_) = 0;

   /// Perform the operation: @a y_ = @a alpha * @a x_ + @a beta * @a y_.
   /// @see braid_PtFcnSum.
   virtual braid_Int Sum(braid_Real   alpha,
                         braid_Vector x_,
                         braid_Real   beta,
                         braid_Vector y_) = 0;

   /// Compute in @a *norm_ptr an appropriate spatial norm of @a u_.
   /// @see braid_PtFcnSpatialNorm.
   virtual braid_Int SpatialNorm(braid_Vector  u_,
                                 braid_Real   *norm_ptr) = 0;

   /// @see braid_PtFcnAccess.
   virtual braid_Int Access(braid_Vector       u_,
                            BraidAccessStatus &astatus) = 0;

   /// @see braid_PtFcnBufSize.
   virtual braid_Int BufSize(braid_Int         *size_ptr,
                             BraidBufferStatus &bstatus) = 0;

   /// @see braid_PtFcnBufPack.
   virtual braid_Int BufPack(braid_Vector       u_,
                             void              *buffer,
                             BraidBufferStatus &bstatus) = 0;

   /// @see braid_PtFcnBufUnpack.
   virtual braid_Int BufUnpack(void              *buffer,
                               braid_Vector      *u_ptr,
                               BraidBufferStatus &bstatus) = 0;

   /// This function may be optionally defined by the user. This provides
   /// a once-per-processor access to the user's app function at certain
   /// points during the code (see documentation for more details).
   /// To turn on sync, use core.SetSync()
   /// @see braid_PtFcnSync.
   virtual braid_Int Sync(BraidSyncStatus &sstatus)
   {
      fprintf(stderr, "Braid C++ Wrapper Warning: turn off sync "
              "until the Sync function been user implemented\n");
      return 0;
   }
};

// Wrapper for BRAID's TriStatus
class BraidTriStatus
{
  private:
      braid_TriStatus pstatus;
  
   public:
      BraidTriStatus(braid_TriStatus _pstatus)
      {
         pstatus = _pstatus;
      }

      void GetTriT(braid_Real *t_ptr, braid_Real *tprev_ptr, braid_Real *tnext_ptr) {
        braid_TriStatusGetTriT(pstatus, t_ptr, tprev_ptr, tnext_ptr);
      }

      void GetIter(braid_Int *iter_ptr) {
        braid_TriStatusGetIter(pstatus, iter_ptr);
      }
        
      void GetLevel(braid_Int *level_ptr) {
        braid_TriStatusGetLevel(pstatus, level_ptr);
      }
      
      void GetNLevels(braid_Int *nlevels_ptr) {
        braid_TriStatusGetNLevels(pstatus, nlevels_ptr);
      }

      void GetNRefine(braid_Int *nrefine_ptr) {
        braid_TriStatusGetNRefine(pstatus, nrefine_ptr);
      }
    
      void GetNTPoints(braid_Int *ntpoints_ptr) {
        braid_TriStatusGetNTPoints(pstatus, ntpoints_ptr);
      }

      void GetOldFineTolx(braid_Real *old_fine_tolx_ptr) {
        braid_TriStatusGetOldFineTolx(pstatus, old_fine_tolx_ptr);
      }
    
      void GetRNorms(braid_Int *nrequest_ptr, braid_Real *rnorms_ptr) {
        braid_TriStatusGetRNorms(pstatus, nrequest_ptr, rnorms_ptr);
      }

      void GetT(braid_Real *t_ptr) {
        braid_TriStatusGetT(pstatus, t_ptr);
      }

      void GetTIndex(braid_Int *index_ptr) {
        braid_TriStatusGetTIndex(pstatus, index_ptr);
      }

      void GetTstartTstop(braid_Real *tstart_ptr, braid_Real *tstop_ptr) {
        braid_TriStatusGetTstartTstop(pstatus, tstart_ptr, tstop_ptr);
      }
    
      void GetTstop(braid_Real *tstop_ptr) {
        braid_TriStatusGetTstop(pstatus, tstop_ptr);
      }

      void GetXRelax(braid_Int *xrelax_ptr) {
        braid_TriStatusGetXRelax(pstatus, xrelax_ptr);
      }

      void GetTol(braid_Real *tol_ptr) {
        braid_TriStatusGetTol(pstatus, tol_ptr);
      }

      void SetOldFineTolx(braid_Real tolx) {
        braid_TriStatusSetOldFineTolx(pstatus, tolx);
      }

      void SetRFactor(braid_Real rfactor) {
        braid_TriStatusSetRFactor(pstatus, rfactor);
      }
      
      void SetRSpace(braid_Real rspace) {
        braid_TriStatusSetRSpace(pstatus, rspace);
      }

      void SetTightFineTolx(braid_Real tft) {
        braid_TriStatusSetTightFineTolx(pstatus, tft);
      }
};

// TriMGRIT specifics
static braid_Int _BraidTriResidual(braid_App     _app,
                                   braid_Vector  _uleft,
                                   braid_Vector  _uright,
                                   braid_Vector  _f,
                                   braid_Vector  _r,
                                   braid_TriStatus _pstatus)
{
  TriBraidApp *app = (TriBraidApp *)_app;
  BraidTriStatus pstatus(_pstatus);
  return app->TriResidual(_uleft, _uright, _f, _r, pstatus);
}


static braid_Int _BraidTriSolve(braid_App      _app,
                                braid_Vector   _uleft,
                                braid_Vector   _uright,
                                braid_Vector   _fleft,
                                braid_Vector   _fright,
                                braid_Vector   _f,
                                braid_Vector   _u,
                                braid_TriStatus _pstatus)
{
  TriBraidApp *app = (TriBraidApp *)_app;
  BraidTriStatus pstatus(_pstatus);
  return app->TriSolve(_uleft, _uright, _fleft, _fright, _f, _u, pstatus);
}


static braid_Int _TriBraidAppClone(braid_App     _app,
                                braid_Vector  _u,
                                braid_Vector *v_ptr)
{
  TriBraidApp *app = (TriBraidApp *)_app;
  return app->Clone(_u, v_ptr);
}


static braid_Int _TriBraidAppInit(braid_App     _app,
                               braid_Real    t,
                               braid_Vector *u_ptr)
{
  TriBraidApp *app = (TriBraidApp *)_app;
  return app->Init(t, u_ptr);
}


static braid_Int _TriBraidAppFree(braid_App    _app,
                               braid_Vector _u)
{
  TriBraidApp *app = (TriBraidApp *)_app;
  return app->Free(_u);
}


static braid_Int _TriBraidAppSum(braid_App    _app,
                              braid_Real   alpha,
                              braid_Vector _x,
                              braid_Real   beta,
                              braid_Vector _y)
{
  TriBraidApp *app = (TriBraidApp *)_app;
  return app->Sum(alpha, _x, beta, _y);
}


static braid_Int _TriBraidAppSpatialNorm(braid_App     _app,
                                      braid_Vector  _u,
                                      braid_Real   *norm_ptr)
{
  TriBraidApp *app = (TriBraidApp *)_app;
  return app->SpatialNorm(_u, norm_ptr);
}


static braid_Int _TriBraidAppAccess(braid_App          _app,
                                 braid_Vector       _u,
                                 braid_AccessStatus _astatus)
{
  TriBraidApp *app = (TriBraidApp *)_app;
  BraidAccessStatus astatus(_astatus);
  return app->Access(_u, astatus);
}

static braid_Int _TriBraidAppBufSize(braid_App  _app,
                                  braid_Int *size_ptr,
                                  braid_BufferStatus _bstatus)
{
  TriBraidApp *app = (TriBraidApp *)_app;
  BraidBufferStatus bstatus(_bstatus);
  return app->BufSize(size_ptr, bstatus);
}


static braid_Int _TriBraidAppBufPack(braid_App     _app,
                                  braid_Vector  _u,
                                  void         *buffer,
                                  braid_BufferStatus  _bstatus)
{
  TriBraidApp *app = (TriBraidApp *)_app;
  BraidBufferStatus bstatus(_bstatus);
  return app->BufPack(_u, buffer, bstatus);
}


static braid_Int _TriBraidAppBufUnpack(braid_App     _app,
                                    void         *buffer,
                                    braid_Vector *u_ptr,
                                    braid_BufferStatus _bstatus)
{
  TriBraidApp *app = (TriBraidApp *)_app;
  BraidBufferStatus bstatus(_bstatus);
  return app->BufUnpack(buffer, u_ptr, bstatus);
}

// Wrapper for BRAID's core object
class BraidTriCore {
private:
   braid_Core core;

public:
  BraidTriCore(MPI_Comm comm_world, TriBraidApp *app) {
    braid_InitTriMGRIT(comm_world, app->comm_t, app->tstart, app->tstop,
                       app->ntime, (braid_App)app, _BraidTriResidual,
                       _BraidTriSolve, _TriBraidAppInit, _TriBraidAppClone,
                       _TriBraidAppFree, _TriBraidAppSum, _TriBraidAppSpatialNorm,
                       _TriBraidAppAccess, _TriBraidAppBufSize, _TriBraidAppBufPack,
                       _TriBraidAppBufUnpack, &core);
  }

   void SetMaxLevels(braid_Int max_levels) { braid_SetMaxLevels(core, max_levels); }

   void SetIncrMaxLevels() { braid_SetIncrMaxLevels(core); }

   void SetSkip(braid_Int skip) { braid_SetSkip(core, skip); }

   void SetMinCoarse(braid_Int min_coarse) { braid_SetMinCoarse(core, min_coarse); }

   void SetNRelax(braid_Int level, braid_Int nrelax)
   { braid_SetNRelax(core, level, nrelax); }

   void SetAbsTol(braid_Real tol) { braid_SetAbsTol(core, tol); }

   void SetRelTol(braid_Real tol) { braid_SetRelTol(core, tol); }

   void SetTemporalNorm(braid_Int tnorm) { braid_SetTemporalNorm(core, tnorm); }

   void SetCFactor(braid_Int level, braid_Int cfactor)
   { braid_SetCFactor(core, level, cfactor); }

   /// If cfactor0 > -1, set the cfactor for level 0 to cfactor0.
   void SetAggCFactor(braid_Int cfactor0)
   {
      if (cfactor0 > -1)
         braid_SetCFactor(core, 0, cfactor0);

    /* Use cfactor0 on all levels until there are < cfactor0 points
       on each processor. */
    //BraidApp *app = (BraidApp *) core->app;
    //braid_Int nt = app->ntime, pt;
    //MPI_Comm_size(app->comm_t, &pt);
    //if (cfactor0 > -1)
    //{
    //   braid_Int level = (braid_Int) (log10((nt + 1) / pt) / log10(cfactor0));
    //   for (braid_Int i = 0; i < level; i++)
    //      braid_SetCFactor(core, i, cfactor0);
    //}
   }

   void SetSpatialCoarsenAndRefine()
   {
      braid_SetSpatialCoarsen(core, _BraidAppCoarsen);
      braid_SetSpatialRefine(core, _BraidAppRefine);
   }

   void SetSync() { braid_SetSync(core, _BraidAppSync); }

   void SetResidual() { braid_SetResidual(core, _BraidAppResidual); }

   void SetMaxIter(braid_Int max_iter) { braid_SetMaxIter(core, max_iter); }

   void SetPrintLevel(braid_Int print_level) { braid_SetPrintLevel(core, print_level); }

   void SetSeqSoln(braid_Int use_seq_soln) { braid_SetSeqSoln(core, use_seq_soln); }

   void SetPrintFile(const char *printfile_name) { braid_SetPrintFile(core, printfile_name); }

   void SetAccessLevel(braid_Int access_level) { braid_SetAccessLevel(core, access_level); }

   void SetFMG() { braid_SetFMG(core); }

   void SetNFMG(braid_Int k) { braid_SetNFMG(core, k); }

   void SetNFMGVcyc(braid_Int nfmg_Vcyc) { braid_SetNFMGVcyc(core, nfmg_Vcyc); }

   void SetStorage(braid_Int storage) { braid_SetStorage(core, storage); }

   void SetRefine(braid_Int refine) {braid_SetRefine(core, refine);}

   void SetMaxRefinements(braid_Int max_refinements) {braid_SetMaxRefinements(core, max_refinements);}

   void GetNumIter(braid_Int *niter_ptr) { braid_GetNumIter(core, niter_ptr); }

   void GetRNorms(braid_Int *nrequest_ptr, braid_Real *rnorms) { braid_GetRNorms(core, nrequest_ptr, rnorms); }
   
   void GetNLevels(braid_Int *nlevels_ptr) { braid_GetNLevels(core, nlevels_ptr); }

   void Drive() { braid_Drive(core); }

   ~BraidTriCore() { braid_Destroy(core); }
};

#endif
