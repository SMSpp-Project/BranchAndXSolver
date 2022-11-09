/*--------------------------------------------------------------------------*/
/*---------------------- File BranchAndXSolver.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the BranchAndXSolver class.
 *
 * \author Federica Di Pasquale \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 * 
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Federica Di Pasquale, Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "BranchAndXSolver.h"

/*--------------------------------------------------------------------------*/
/*--------------------------------- MACROS ---------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE AND USING --------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------- CONSTANTS --------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register BranchAndXSolver to the Solver factory

SMSpp_insert_in_factory_cpp_0( BranchAndXSolver );

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS OF BranchAndXSolver ------------------------*/
/*--------------------------------------------------------------------------*/

void BranchAndXSolver::set_Block( Block * block ) {

 if( block == f_Block )       // nothing to do        
  return;
 
 Solver::set_Block( block );  // attach to the new Block

 f_Solver->set_Block( f_Block );
 f_sense = f_Block->get_objective_sense(); 
 }

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
/// 

int BranchAndXSolver::compute( bool changedvars ) {
 
 lock();                    // lock the mutex

 process_outstanding_Modification();
 
 /*-------------------------------------------------------------------------*/
 
 // Start Timer
 auto start = std::chrono::system_clock::now();

 // The root node is the original problem
 Node * root = new Node( f_Block->get_R3_Block() ); 
 
 // Queue of subproblems still to be processed
 Queue Q( root );

 // Number of iterations
 Index n_iter = 0;

 while( ! Q.empty() ) {
  
  // Select next subproblem
  auto node = Q.pop();

  // Compute the relaxation 
  f_Solver->set_Block( node->block() );
  f_Solver->compute();

  // Update lower and upper bounds - - - - - - - - - - - - - - - - - - - - - -
  
  if( get_objective_sense() == Objective::eMax ) {      // MAX case
   node->set_ub( f_Solver->get_var_value() );           // upper bound
   node->set_lb( f_RelaxationSolver->get_true_lb() );   // lower bound
   f_lb = std::max( node->get_lb() , f_lb );            // global lower bound
   }

  if( get_objective_sense() == Objective::eMin ) {      // MIN case
   node->set_lb( f_Solver->get_var_value() );           // lower bound
   node->set_ub( f_RelaxationSolver->get_true_ub() );   // upper bound
   f_ub = std::min( node->get_ub() , f_ub );            // global upper bound
   }

  // Pruning Rules - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // Check feasibility
  if( node->is_feasible( RelAcc ) ) {
   delete node;
   continue;
   } 

  // Pruning by bound
  if( get_objective_sense() == Objective::eMax && node->get_ub() < f_lb ) {
   delete node;
   continue;
   } 
  if( get_objective_sense() == Objective::eMin && node->get_lb() > f_ub ) {
   delete node;
   continue;
   } 

  // Branch - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  auto branches = f_RelaxationSolver->branch();
  
  for( auto branch : branches ) {
   auto block = node->block()->get_R3_Block();
   branch->apply( block ); 
   Q.push( new Node( block ) );
   f_n_nodes++;
   delete branch;                  
   }

   // Check time/iteration limits - - - - - - - - - - - - - - - - - - - - - -

   if( n_iter > MaxIter )
    break;
   n_iter++;

   auto end = std::chrono::system_clock::now();
   std::chrono::duration< double > elapsed = end - start;

   if( elapsed.count() > MaxTime )
    break;

 }

 /*-------------------------------------------------------------------------*/ 

 f_obj = ( get_objective_sense() == Objective::eMax ) ? f_lb : f_ub;

 Return_OK:
 
 unlock();                  // unlock the mutex     
 
 return( kOK );
 
 } // end( BranchAndXSolver::compute() )

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------- End File BranchAndXSolver.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
