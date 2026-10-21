/*--------------------------------------------------------------------------*/
/*----------------------- File BranchAndXSolver.cppm -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Interface of the C++20 module smspp.branchandxsolver, which exports every
 * name that the headers of BranchAndXSolver declare at namespace scope,
 * together with the modules it depends on; a translation unit can then write
 *
 *     import smspp.branchandxsolver;
 *
 * in place of the #include of the headers. The macros (e.g., those of the
 * factory in SMSFactory.h) are not exported by a module, so a file that uses
 * them still includes the header that defines them.
 *
 * This file is written by SMS++/cmake/smspp_module_exports.py, which reads
 * the names from the AST of the headers: it is not to be edited by hand.
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
 * \copyright &copy; by Antonio Frangioni, Filippo Magi, Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*------------------------- GLOBAL MODULE FRAGMENT -------------------------*/
/*--------------------------------------------------------------------------*/

module;

#include "BranchAndXSolver.h"
#include "BranchingRule.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- MODULE PURVIEW -----------------------------*/
/*--------------------------------------------------------------------------*/

export module smspp.branchandxsolver;

export import smspp;

export namespace SMSpp_di_unipi_it {

 // BranchAndXSolver.h
 using SMSpp_di_unipi_it::BranchAndXSolver;
 using SMSpp_di_unipi_it::DFSNode;
 using SMSpp_di_unipi_it::ExploringNode;
 using SMSpp_di_unipi_it::Node;
 using SMSpp_di_unipi_it::OpenList;

 // BranchingRule.h
 using SMSpp_di_unipi_it::BranchCandidate;
 using SMSpp_di_unipi_it::BranchingFeatures;
 using SMSpp_di_unipi_it::LearnedBranching;
 using SMSpp_di_unipi_it::StrongBranching;
}

/*--------------------------------------------------------------------------*/
/*--------------------- End File BranchAndXSolver.cppm ---------------------*/
/*--------------------------------------------------------------------------*/
