/*--------------------------------------------------------------------------*/
/*------------------------- File ParallelSolver.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *abstract* template class ParallelSolver, a Solver
 * able to work in parallel with clones of itself: each ParallelSolver can
 * create clones via clone(); the clones are not attached to the Block, but
 * the master solver (the one that is) propagates to them every Modification
 * it receives. This is meant to allow multiple solvers to work in parallel
 * on the same problem, e.g., inside a parallel Branch-and-Bound. WARNING:
 * WORK IN PROGRESS - not used by BranchAndXSolver yet.
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
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __ParallelSolver
 #define __ParallelSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>
#include <list>
#include <memory>

#if defined( __cpp_lib_parallel_algorithm )
 #include <execution>
#endif

#include "Solver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS ParallelSolver ----------------------------*/
/*--------------------------------------------------------------------------*/
/// a Solver able to work in parallel with clones of itself
/** The template class ParallelSolver (a CRTP base: T must derive from
 * ParallelSolver< T >) provides the machinery to create clones of a Solver
 * that share the master's view of the Block: the clones are not attached to
 * the Block, the master propagates to them every Modification it receives.
 */

template< typename T >
class ParallelSolver : public Solver {

 static_assert( std::is_base_of_v< ParallelSolver< T > , T > ,
                "ParallelSolver: T must be derived from ParallelSolver< T >"
                );

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

 ParallelSolver() : Solver() , master( nullptr ) {}    ///< constructor

 /// destructor: detaches the clone from its master

 ~ParallelSolver() override {
  if( master )
   master->removeSon( this );
  }

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED METHODS ----------------------------*/
/*--------------------------------------------------------------------------*/

 /// create an independent copy of this ParallelSolver
 /** Create an independent copy of this ParallelSolver. The clone is
  * registered as a son of the master of this solver (or of this solver
  * itself, if it is the master): it is not attached to the Block, and it
  * receives from its master all the Modification added to it. */

 std::unique_ptr< T > clone( void ) {
  auto newSolver = std::make_unique< T >( static_cast< const T & >( *this )
                                          );
  if( master ) {
   master->addSon( newSolver.get() );
   newSolver->master = master;
   }
  else {
   addSon( newSolver.get() );
   newSolver->master = static_cast< T * >( this );
   }
  return( newSolver );
  }

/*--------------------------------------------------------------------------*/
 /// notify this solver that a new son solver has been created

 void addSon( ParallelSolver * son ) {
  if( master )
   master->sonSolvers.push_back( son );
  else
   sonSolvers.push_back( son );
  son->master = this;
  }

/*--------------------------------------------------------------------------*/
 /// notify this solver that a son solver has been deleted

 void removeSon( ParallelSolver * son ) {
  if( master )
   master->sonSolvers.remove( son );
  else
   sonSolvers.remove( son );
  son->master = nullptr;
  }

/*--------------------------------------------------------------------------*/
 /// receive a Modification and propagate it to the son solvers

 void add_Modification( sp_Mod & mod ) override {
  Solver::add_Modification( mod );
  if( ! master )            // only the master propagates the Modification
   parallel_for_each( sonSolvers.begin() , sonSolvers.end() ,
                      [ & ]( Solver * s ) { s->add_Modification( mod ); } );
  }

/*--------------------------------------------------------------------------*/
 /// std::for_each, parallel whenever the standard library allows it

 template< typename Iterator , typename Func >
 void parallel_for_each( Iterator first , Iterator last , Func f ) {
  #if defined( __cpp_lib_parallel_algorithm )
   std::for_each( std::execution::par_unseq , first , last , f );
  #else
   std::for_each( first , last , f );
  #endif
  }

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 T * master;                  ///< the master solver, attached to the Block

 /// the son solvers, fed with the Modification by the master
 std::list< ParallelSolver * > sonSolvers;

/*--------------------------------------------------------------------------*/

 };  // end( class( ParallelSolver ) )

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* ParallelSolver.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File ParallelSolver.h -------------------------*/
/*--------------------------------------------------------------------------*/
