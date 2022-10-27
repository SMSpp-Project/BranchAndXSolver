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

 f_RelaxationSolver->set_Block( f_Block );
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
 
 // The root node is the original problem
 Node * root = new Node( f_Block->get_R3_Block() ); 
 
 // Queue of subproblems still to be processed
 Queue Q( root );

 while( ! Q.empty() ) {
  
  // Select next subproblem
  auto node = Q.pop();

  // Compute the relaxation and update node bound (lb or ub)
  f_RelaxationSolver->set_Block( node->block() );
  f_RelaxationSolver->compute();
  node->set_bound( f_RelaxationSolver->get_var_value() ); 

  // Pruning Rules- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  double gap = + Inf< double >();  // to check if the problem has been solved

  if( f_sense == Objective::eMax ) {    // MAX                           
   if( RlxSlv->get_true_lb() > f_lb )
    f_lb = RlxSlv->get_true_lb(); 
   gap = std::abs( RlxSlv->get_true_lb() - node->bound() );        
   }
  else {                                // MIN                                   
   if( RlxSlv->get_true_ub() < f_ub )
    f_ub = RlxSlv->get_true_ub();
   gap = std::abs( RlxSlv->get_true_ub() - node->bound() ); 
   }
  
  // Check feasibility
  if( gap < 1e-04 ) {
   delete node;
   continue;
   }

  // Pruning by bound
  if( ( f_sense == Objective::eMax ) && ( node->bound() < f_lb ) ) {     
   delete node; 
   continue;
   }
  if( ( f_sense == Objective::eMin ) && ( node->bound() > f_ub ) ) {    
   delete node; 
   continue;
   }
  
  // Branch - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  std::vector< Change * > children = RlxSlv->branch();
  
  for( auto child : children ) {
   auto block = node->block()->get_R3_Block();
   child->apply( block ); 
   Q.push( new Node( block ) );
   f_n_nodes++;                  
   }

 }

 /*-------------------------------------------------------------------------*/ 

 f_obj = ( f_sense == Objective::eMax ) ? f_lb : f_ub;

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
