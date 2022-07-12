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
 
 // Queue of subproblems still to be processed
 std::priority_queue< BaXnode * > Q;  

 // The root node is the original problem
 BaXnode * root = new BaXnode(); 

 // Initialize the queue with the original problem   
 Q.push( root );
 
 // last processed node
 BaXnode * old_node = nullptr;

 // for each subproblem in the queue Q
 for( ; ! Q.empty() ; Q.pop() ) {

  // select new subproblem
  auto node = Q.top();

  // update current subproblem 
  update_subproblem( old_node , node );
  
  // compute relaxation
  f_RelaxationSolver->compute();

  // update bounds of node
  node->set_bound( f_RelaxationSolver->get_var_value() ); // lb o ub

  // Pruning Rules- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  double gap = + Inf< double >();  // to check if the problem has been solved

  if( f_sense == Objective::eMin ) {       // if it is a minimization problem
   f_ub = std::min( f_ub , RlxSlv->get_true_ub() );     // global upper bound
   gap = std::abs( RlxSlv->get_true_ub() - node->bound() );        // gap
   if( node->bound() >= f_ub ){
    //std::cout << "pruned by bound\n";
    continue; 
    }
   }
  else
  if( f_sense == Objective::eMax ) {       // if it is a maximization problem
   f_lb = std::max( f_lb , RlxSlv->get_true_lb() );     // global lower bound 
   gap = std::abs( RlxSlv->get_true_lb() - node->bound() );        // gap
   if( node->bound() <= f_lb ){
    //std::cout << "pruned by bound\n";
    continue; 
    }
   }
   
   // check if the problem has been solved to optimality
   if( gap < 1e-04 ) {
    //std::cout << "pruned by feasibility\n";
    continue;   
    }


  // Branch - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  std::vector< Change * > children = RlxSlv->branch();
  for( auto child : children ) {
   auto new_node = new BaXnode( node , node->depth() + 1 , child );  
   Q.push( new_node );                  
   node->add_child( new_node );
   }

  // Stopping conditions
  // accuracy | max_iter | ...
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

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// update current subproblem and its path and reverse path from the root

 void BranchAndXSolver::update_subproblem( BaXnode * old_node , 
                                           BaXnode * node ) {

  // if it is the root node
  if( ! node->father() )
   return;

  // compute least common ancestor (lca)
  std::vector< Change * > path;

  // apply undoChg up to lca

  // apply chg from lca to the current node - 1

  // store last undoChg with node->set_undoChg()
 
  old_node = node;
  }
/*--------------------------------------------------------------------------*/
/*--------------------- End File BranchAndXSolver.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
