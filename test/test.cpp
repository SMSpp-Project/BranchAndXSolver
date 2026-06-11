/*--------------------------------------------------------------------------*/
/*----------------------------- File test.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Minimal driver for BranchAndXSolver: solves a tiny mixed Binary
 * Knapsack instance with the GreedyRelaxationSolver continuous relaxation
 * under each of the three tree exploration strategies (DFS, BFS, BestFS)
 * and prints the optimal values.
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
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>

#include "BranchAndXSolver.h"

#include "GreedyRelaxationSolver.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- USING ------------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- MAIN ------------------------------------*/
/*--------------------------------------------------------------------------*/

int main()
{
 auto block = new BinaryKnapsackBlock();
 std::vector< double > weights = { 1.5 , 2 , 3 , 4 , 5.5 };
 std::vector< double > profits = { 10.1 , 21.0 , 32.0 , 43.0 , 51.0 };
 block->load( 5 , 8 , std::move( weights ) , std::move( profits ) );

 auto solver = new GreedyRelaxationSolver();
 BranchAndXSolver tree( solver , nullptr , block );

 tree.set_par( BranchAndXSolver::intSolveMethod ,
               BranchAndXSolver::DFS );
 std::cout << "Solving DFS:" << std::endl;
 std::cout << "Return value: " << tree.compute() << std::endl;
 std::cout << "Best solution found: " << tree.get_var_value() << std::endl;

 tree.set_par( BranchAndXSolver::intSolveMethod ,
               BranchAndXSolver::BFS );
 std::cout << "Solving BFS:" << std::endl;
 std::cout << "Return value: " << tree.compute() << std::endl;
 std::cout << "Best solution found: " << tree.get_var_value() << std::endl;

 tree.set_par( BranchAndXSolver::intSolveMethod ,
               BranchAndXSolver::BestFS );
 std::cout << "Solving Best First:" << std::endl;
 std::cout << "Return value: " << tree.compute() << std::endl;
 std::cout << "Best solution found: " << tree.get_var_value() << std::endl;

 return( 0 );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------- End File test.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
