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

 // an (un)fixing Change is applied to the internal mirror only, reading the
 // acted-upon items via the uniform polymorphic view of the Change (see
 // BinaryKnapsackBlockChange::num_items() / item()): the concrete (ranged /
 // subset) type of the Change is never needed
 const bool fixing = ( CHG->type() == BinaryKnapsackBlockChange::eFixX );
 const Index n = CHG->num_items();
 const auto & data = CHG->data();

 Change * undo = nullptr;
 if( doUndo ) {
  // the undo under the LIFO discipline of BranchAndXSolver: un-fixing the
  // same items, or re-fixing them to the values captured at call time
  std::vector< double > old_data;
  Block::Subset nms;
  old_data.reserve( n );
  nms.reserve( n );
  for( Index k = 0 ; k < n ; ++k ) {
   const Index i = CHG->item( k );
   nms.push_back( i );
   old_data.push_back( v_fxd[ i ] == 2 ? 1 : 0 );
   }
  undo = new BinaryKnapsackBlockSbstChange(
                     fixing ? BinaryKnapsackBlockChange::eUnfixX
                            : BinaryKnapsackBlockChange::eFixX ,
                     std::move( old_data ) , std::move( nms ) );
  }

 for( Index k = 0 ; k < n ; ++k )
  v_fxd[ CHG->item( k ) ] = fixing ? ( data[ k ] == 1 ? 2 : 1 ) : 0;

 return( undo );

 }  // end( GreedyRelaxationSolver::apply )

/*--------------------------------------------------------------------------*/
/*------------------ End File GreedyRelaxationSolver.cpp -------------------*/
/*--------------------------------------------------------------------------*/
