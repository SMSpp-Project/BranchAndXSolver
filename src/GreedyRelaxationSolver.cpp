/*--------------------------------------------------------------------------*/
/*--------------------- File GreedyRelaxationSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the *concrete* class GreedyRelaxationSolver, the
 * RelaxationSolver [see ChangeSolver.h] for the continuous relaxation of
 * the Binary Knapsack problem in view of its use within
 * BranchAndXSolver: everything except apply() is inherited from
 * GreedyRelaxationBinaryKnapsackSolver of the BinaryKnapsackBlock module.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Filippo Magi \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni, Filippo Magi, Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "GreedyRelaxationSolver.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using Index = Block::Index;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register GreedyRelaxationSolver to the Solver factory

SMSpp_insert_in_factory_cpp_1( GreedyRelaxationSolver );

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

Change * GreedyRelaxationSolver::apply( Change * chg , bool doUndo )
{
 auto CHG = dynamic_cast< BinaryKnapsackBlockChange * >( chg );

 if( ! CHG )
  throw( std::invalid_argument( "GreedyRelaxationSolver::apply: the Change "
                                "must be a BinaryKnapsackBlockChange" ) );

 if( ( CHG->type() != BinaryKnapsackBlockChange::eFixX ) &&
     ( CHG->type() != BinaryKnapsackBlockChange::eUnfixX ) )
  return( CHG->apply( f_Block , doUndo ) );   // not an (un)fix: to the Block

 // an (un)fixing Change is applied to the internal mirror only [see
 // BinaryKnapsackSolver::apply_to_mirror()]: the next compute() re-solves
 // the relaxation under the new fixings without any access to the Block
 return( apply_to_mirror( CHG , doUndo ) );

 }  // end( GreedyRelaxationSolver::apply )

/*--------------------------------------------------------------------------*/
/*------------------ End File GreedyRelaxationSolver.cpp -------------------*/
/*--------------------------------------------------------------------------*/
