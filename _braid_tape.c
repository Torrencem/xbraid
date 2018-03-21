/** \file _braid_tape.c
 * \brief Source code for tape routines. See _braid_tape.h for more documentation.
 *
 */

#include "_braid.h"
#include "_braid_tape.h"
#include "_braid_base.h"

#ifndef DEBUG
#define DEBUG 0
#endif


void 
_braid_TapeInit(_braid_Tape* head)
{
   head = NULL;
}

_braid_Tape* 
_braid_TapePush(_braid_Tape* head, void* data_ptr)
{
   _braid_Tape* tmp = (_braid_Tape*)malloc(sizeof(_braid_Tape));
   if (tmp == NULL)
   {
      printf("MALLOC ERROR!\n") ;
      exit(1);
   }
   if ( _braid_TapeIsEmpty(head) ) tmp->size = 1;
   else tmp->size = head->size + 1;
   tmp->data_ptr = data_ptr;
   tmp->next = head;
   head = tmp;

   return head;
}

_braid_Tape* 
_braid_TapePop(_braid_Tape *head)
{
    _braid_Tape* tmp = head;
    head->size = tmp->size--;
    head = head->next;
    free(tmp);

    return head;
}

braid_Int 
_braid_TapeIsEmpty(_braid_Tape* head)
{
    return head == NULL ? 1 : 0;
}

braid_Int
_braid_TapeGetSize(_braid_Tape* head)
{
    if ( _braid_TapeIsEmpty(head) ) return 0;
    else return head->size;
}


void 
_braid_TapeDisplayBackwards(braid_Core core, _braid_Tape* head, void (*displayfct)(braid_Core core, void* data_ptr))
{
    _braid_Tape* current;
    current = head;
    if (current!=NULL)
    {
        do 
        {
            /* Call the display function */
            (*displayfct)(core, current->data_ptr);
            /* Move to next element */
            current = current->next;
        }
        while (current!=NULL);
    }
    else
    {
        printf("Tape is empty\n");
    }
  
}


void 
_braid_TapeEvaluate(braid_Core core)
{
    /* Get the tape */
    _braid_Tape* actiontape = _braid_CoreElt(core, actiontape);

    while ( !_braid_TapeIsEmpty(actiontape) )
    {
       /* Get the action */
       _braid_Action* action = (_braid_Action*) actiontape->data_ptr;

       /* Call the differentiated action */
       _braid_DiffCall(action);

       /* Pop the action from the tape */
       actiontape = _braid_TapePop( actiontape );

       /* Free memory of the action */
       free(action);
    }
    
    /* Update the actiontape in the core */
    _braid_CoreElt(core, actiontape) = actiontape;
}

void
_braid_DiffCall(_braid_Action* action)
{
   
   /* Call the corresponding differentiated action */
   switch (action->braidCall)
   {
      case STEP : 
      {
         _braid_BaseStep_diff(action);
         break;
      }
      case INIT: 
      {
         /* TODO: so something, if design is part of  initial conditions!!  */
         break;
      }
      case CLONE: 
      {
         _braid_BaseClone_diff(action);
         break;
      }
      case FREE: 
      {
         /* Do nothing */
         break;
      }
      case SUM: 
      {
         _braid_BaseSum_diff(action);
         break;
      }
      case ACCESS: 
      {
         /* Do nothing (because access is only for output!) */
         break;
      }
      case BUFPACK: 
      {
         _braid_BaseBufPack_diff(action, _braid_CoreElt(action->core, app));
         break;
      }
      case BUFUNPACK: 
      {
         _braid_BaseBufUnpack_diff(action, _braid_CoreElt(action->core, app));
         break;
      } 
      case OBJECTIVET:
      {
          /* TODO: Call the differentiated routine */
          _braid_BaseObjectiveT_diff(action);
          break;
      }
   }
}


void
_braid_TapeSetSeed(braid_Core core)
{
   _braid_Grid *fine_grid;
   braid_BaseVector u_out;
   braid_Int clower, cfactor, iclocal, iupper, sflag;
   braid_Optim optim;
   braid_App app;

   fine_grid = _braid_CoreElt(core, grids)[0];
   clower    = _braid_GridElt(fine_grid, clower);
   iupper    = _braid_GridElt(fine_grid, iupper);
   cfactor   = _braid_GridElt(fine_grid, cfactor);
   optim     = _braid_CoreElt(core, optim);
   app       = _braid_CoreElt(core, app);

   /* Loop over all coarse points on finest grid level */
   for (braid_Int ic=clower; ic <= iupper; ic += cfactor)
   {
      /* Get the output vector u and its index in the ua vector  */
      _braid_UGetIndex(core, 0, ic, &iclocal, &sflag);
      _braid_UGetVectorRef(core, 0, ic, &u_out);

      /* Set the seed using the optimization adjoints */
      _braid_CoreFcn(core, sum)(app, 1.0, optim->adjoints[iclocal]->userVector, 0.0, u_out->adjoint->userVector);

      /* Debug: */
      // printf("Seeding for %d %p\n", iclocal, u_out->adjoint );
      // _braid_CoreFcn(core, access)(app, u_out->adjoint->userVector, NULL);

   }
}

void 
_braid_TapeResetInput(braid_Core core)
{
   braid_BaseVector u_out;
   braid_Adjoint    adjoint_copy;
   _braid_Grid     *fine_grid;
   braid_Int        clower, iupper, cfactor, iclocal, sflag, ic;
   // braid_Optim      optim;

   fine_grid = _braid_CoreElt(core, grids)[0];
   clower    = _braid_GridElt(fine_grid, clower);
   iupper    = _braid_GridElt(fine_grid, iupper);
   cfactor   = _braid_GridElt(fine_grid, cfactor);
   // optim     = _braid_CoreElt(core, optim);
  
   /* Loop over all C-points */
   for (ic=clower; ic <= iupper; ic += cfactor)
   {
      /* Get the vector in ua and its local index */
      _braid_UGetIndex(core, 0, ic, &iclocal, &sflag);
      _braid_UGetVectorRef(core, 0, ic, &u_out);

      /* Copy a pointer to its adjoint to the tapeinput vector */
      _braid_AdjointCopy(u_out->adjoint, &adjoint_copy );
      _braid_CoreElt(core, optim)->tapeinput[iclocal] = adjoint_copy;
   }
}

// void
// _braid_TapeDisplayAction(braid_Core core,void* data_ptr){

//     /* Get the action */    
//     _braid_Action* action = (_braid_Action*)(data_ptr);
//     /* Print action information */
//     printf("%d: %s ", action->myid, _braid_CallGetName(action->braidCall));
//     printf("\n");

// }

// void
// _braid_TapeDisplayPrimal(braid_Core core,void* data_ptr)
// {
//     /* Get the braid_vector */
//    braid_Vector vector = (braid_Vector) (data_ptr);

//    /*--- Display the vector --*/
//    braid_AccessStatus   astatus = (braid_AccessStatus)core;
//    _braid_CoreFcn(core, access)(_braid_CoreElt(core, app), vector->primal, astatus );
// }

// void
// _braid_TapeDisplayInt(braid_Core core,void* data_ptr)
// {
//     /* Get the integer*/
//     int* int_ptr = (int* ) data_ptr;

//    /*--- Display the integer --*/
//     printf(" Integer: %d", (*int_ptr) );
// }


const char* _braid_CallGetName(_braid_Call call)
{
    switch (call)
    {
        case STEP:       return "STEP";
        case INIT:       return "INIT";
        case CLONE:      return "CLONE";
        case FREE:       return "FREE";
        case SUM:        return "SUM";
        case BUFPACK:    return "BUFPACK";
        case BUFUNPACK:  return "BUFUNPACK";
        case ACCESS:     return "ACCESS";
        case OBJECTIVET: return "OBJECTIVET";
    }
    return "Not a _braid_Call!";
}