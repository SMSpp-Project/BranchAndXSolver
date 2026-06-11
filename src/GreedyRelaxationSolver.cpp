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
  return( CHG->apply( f_Block , doUndo ) );

 const bool fixing = ( CHG->type() == BinaryKnapsackBlockChange::eFixX );

 // the new value of v_fxd for an index being (un)fixed: 2 = fixed to 1,
 // 1 = fixed to 0, 0 = free
 auto new_fxd = [ & ]( double value ) -> unsigned char {
  return( fixing ? ( value == 1 ? 2 : 1 ) : 0 );
  };

 // the (captured) data for the undo Change of an index: the value it is
 // currently fixed to when un-doing an unfix, immaterial when un-doing a fix
 auto old_val = [ & ]( Index i ) -> double {
  return( v_fxd[ i ] == 2 ? 1 : 0 );
  };

 // ranged change (the branching case) - - - - - - - - - - - - - - - - - - - -
 if( auto change = dynamic_cast< BinaryKnapsackBlockRngdChange * >( CHG ) ) {
  Block::Range rng = change->rng();
  const auto & data = change->data();

  Change * undo = nullptr;
  if( doUndo ) {
   std::vector< double > old_data;
   old_data.reserve( rng.second - rng.first );
   for( Index i = rng.first ; i < rng.second ; ++i )
    old_data.push_back( old_val( i ) );
   undo = new BinaryKnapsackBlockRngdChange(
                     fixing ? BinaryKnapsackBlockChange::eUnfixX
                            : BinaryKnapsackBlockChange::eFixX ,
                     std::move( old_data ) , Block::Range( rng ) );
   }

  for( Index i = rng.first ; i < rng.second ; ++i )
   v_fxd[ i ] = new_fxd( fixing ? data[ i - rng.first ] : 0 );

  return( undo );
  }

 // subset change - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( auto change = dynamic_cast< BinaryKnapsackBlockSbstChange * >( CHG ) ) {
  const auto & subset = change->nms();
  const auto & data = change->data();

  Change * undo = nullptr;
  if( doUndo ) {
   std::vector< double > old_data;
   old_data.reserve( subset.size() );
   for( auto i : subset )
    old_data.push_back( old_val( i ) );
   undo = new BinaryKnapsackBlockSbstChange(
                     fixing ? BinaryKnapsackBlockChange::eUnfixX
                            : BinaryKnapsackBlockChange::eFixX ,
                     std::move( old_data ) , Block::Subset( subset ) );
   }

  Index idx = 0;
  for( auto i : subset )
   v_fxd[ i ] = new_fxd( fixing ? data[ idx++ ] : 0 );

  return( undo );
  }

 throw( std::invalid_argument( "GreedyRelaxationSolver::apply: the Change "
        "must be a BinaryKnapsackBlock[Rngd/Sbst]Change" ) );

 }  // end( GreedyRelaxationSolver::apply )

/*--------------------------------------------------------------------------*/
/*------------------ End File GreedyRelaxationSolver.cpp -------------------*/
/*--------------------------------------------------------------------------*/
