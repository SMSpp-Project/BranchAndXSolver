/*--------------------------------------------------------------------------*/
/*---------------------- File GreedyRelaxationSolver.h ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class GreedyRelaxationSolver, which
 * implements the RelaxationSolver concept [see ChangeSolver.h] for the
 * continuous relaxation of the Binary Knapsack problem represented by a
 * BinaryKnapsackBlock, in view of its use within BranchAndXSolver.
 *
 * All the knapsack machinery - the raw mirror of the instance, the
 * normalized core, the greedy (Dantzig) fractional relaxation, the true
 * lower / upper bounds and the branching on the critical item - is
 * inherited from GreedyRelaxationBinaryKnapsackSolver of the
 * BinaryKnapsackBlock module; this class only adds what the enumeration
 * needs: apply()-ing the (un)fixing Changes produced by branch() *to the
 * internal mirror of the instance* (i.e., without touching the Block, so
 * that the same Solver object can be cheaply moved between the nodes of the
 * tree) and producing the true (rounded) solutions as Solution objects.
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

#ifndef __GreedyRelaxationSolver
 #define __GreedyRelaxationSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "GreedyRelaxationBinaryKnapsackSolver.h"

#include "ChangeSolver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------- CLASS GreedyRelaxationSolver ------------------------*/
/*--------------------------------------------------------------------------*/
/// the RelaxationSolver for the Binary Knapsack continuous relaxation
/** GreedyRelaxationSolver joins GreedyRelaxationBinaryKnapsackSolver (the
 * whole knapsack machinery: instance mirror, normalized core, Dantzig
 * relaxation, true bounds, branching) with the RelaxationSolver interface
 * required by BranchAndXSolver; both hierarchies share the single
 * virtual Solver sub-object. The (un)fixing Changes produced by branch()
 * are applied *internally* [see apply()]: the next compute() re-solves the
 * relaxation under the new fixings without any access to the Block. */

class GreedyRelaxationSolver : public GreedyRelaxationBinaryKnapsackSolver ,
                               public RelaxationSolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

 /// constructor

 GreedyRelaxationSolver() : GreedyRelaxationBinaryKnapsackSolver() ,
                            RelaxationSolver() {}

 ~GreedyRelaxationSolver() override = default;   ///< destructor

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

 /// branch on the critical item (see the BinaryKnapsackBlock module class)

 std::vector< Change * > branch() override {
  return( GreedyRelaxationBinaryKnapsackSolver::branch() );
  }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

 /// valid lower bound on the optimal value of the TRUE problem

 OFValue get_true_lb( void ) override {
  return( GreedyRelaxationBinaryKnapsackSolver::get_true_lb() );
  }

 /// valid upper bound on the optimal value of the TRUE problem

 OFValue get_true_ub( void ) override {
  return( GreedyRelaxationBinaryKnapsackSolver::get_true_ub() );
  }

/*--------------------------------------------------------------------------*/
 /// a true (rounded greedy) solution is always available after compute()

 bool has_true_var_solution( void ) override {
  return( GreedyRelaxationBinaryKnapsackSolver::has_true_var_solution() );
  }

/*--------------------------------------------------------------------------*/
 /// write the rounded greedy solution in the Variable of the Block

 void get_true_var_solution( Configuration * solc = nullptr ) override {
  GreedyRelaxationBinaryKnapsackSolver::get_true_var_solution( solc );
  }

/*--------------------------------------------------------------------------*/
 /// physically construct the true Solution, bypassing the Block

 Solution * get_true_solution( Configuration * solc = nullptr ) override {
  return( new BinaryKnapsackSolution( rounded_x() ) );
  }

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

 /// apply the Change; (un)fixing Changes are applied INTERNALLY
 /** The (un)fixing Changes produced by branch() - and, in general, any
  * BinaryKnapsackBlock[Rngd/Sbst]Change of eFixX / eUnfixX type - are
  * applied to the internal mirror of the instance, NOT to the Block: the
  * next compute() then re-solves the relaxation under the new fixings.
  * Any other Change is forwarded to the Block.
  *
  * With \p doUndo true the returned Change un-does \p chg under the LIFO
  * (depth-first / climb-the-tree) discipline of BranchAndXSolver: the
  * undo of fixing a free variable is unfixing it, and vice-versa with the
  * fixed values captured at call time. */

 Change * apply( Change * chg , bool doUndo = false ) override;

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 SMSpp_insert_in_factory_h;  // insert GreedyRelaxationSolver in the factory

/*--------------------------------------------------------------------------*/

 };  // end( class( GreedyRelaxationSolver ) )

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* GreedyRelaxationSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------------ End File GreedyRelaxationSolver.h ---------------------*/
/*--------------------------------------------------------------------------*/
