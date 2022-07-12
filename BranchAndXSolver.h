/*--------------------------------------------------------------------------*/
/*--------------------- File BranchAndXSolver.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the BranchAndXSolver class, which implements the Solver
 * interface ...
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
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __BranchAndXSolver
 #define __BranchAndXSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <queue>

#include "Block.h"
#include "Change.h"
#include "Objective.h"
#include "RelaxationSolver.h"
#include "Solution.h"
#include "Solver.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup BranchAndXSolver_CLASSES Classes in BranchAndXSolver.h
 *  @{ */


/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS BaXnode --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// 
/**  */

class BaXnode {
 
 public:
 
 using Index = Block::Index;

  // - - - - - - - - - - Constructor and Destructor - - - - - - - - - - - - -
 
  BaXnode( BaXnode * f = nullptr , Index depth = 0 , 
           Change * c = nullptr ) : 
           f_father( f ) , f_depth( depth ) , f_chg( c ) , 
           f_undoChg( nullptr ) ,
           f_bound( + Inf< double >() ) {}
  
  ~BaXnode() = default;
 
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  void set_chg( Change * chg ) { f_chg = chg; }
  
  void set_undoChg( Change * undoChg ) {f_undoChg = undoChg;}
  
  void set_father( BaXnode * father ) { f_father = father; }
  
  void add_child( BaXnode * child ) { v_children.push_back( child ); }
  
  void set_bound( double bound ) { f_bound = bound; }

  void set_depth( Index depth ) { f_depth = depth; }

  Change * chg() { return f_chg; }
  
  Change * undoChg() { return f_undoChg; }
  
  BaXnode * father() { return f_father; }
  
  double bound() { return f_bound; }

  Index depth() { return f_depth; }

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 protected:
 
  Change * f_chg;      // change to apply w.r.t. father node
  Change * f_undoChg;  // undo-change
  BaXnode * f_father;                   // pointer to the father node
  std::vector< BaXnode * > v_children;
  double f_bound;
  Index f_depth;

}; // end( class BaXnode )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS BranchAndXSolver ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// 
/**  */

class BranchAndXSolver : public Solver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *
 *
 *  @{ */

 using Index = Block::Index;
 using Subset = Block::Subset;

/** @} ---------------------------------------------------------------------*/
/*--------------- CONSTRUCTING AND DESTRUCTING BranchAndXSolver ------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing BranchAndXSolver
 *  @{ */

 /// constructor

 BranchAndXSolver() : f_RelaxationSolver( nullptr ) ,
                      f_incumbent( nullptr ) , 
                      f_lb( + Inf< double >() ) ,
                      f_ub( - Inf< double >() ) ,
                      f_obj( + Inf< double >() ) , 
                      f_sense( Objective::eMin ) {}

/*--------------------------------------------------------------------------*/
 /// destructor

 virtual ~BranchAndXSolver() = default;

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *
 *  @{ */

/*--------------------------------------------------------------------------*/
 /// set the (pointer to the) Block that the Solver has to solve
 /** */

 void set_Block( Block * block ) override;

/*--------------------------------------------------------------------------*/
 /// set the (pointer to the) the Relaxation Solver 
 /** */

 void set_RelaxationSolver( Solver * rlxslv ) {
  
  if( f_RelaxationSolver == rlxslv )
   return;
  
  RlxSlv = dynamic_cast< RelaxationSolver * >( rlxslv );
  if( ! RlxSlv )
   throw( std::invalid_argument( "Not able to cast to RelaxationSolver *" ) );

  f_RelaxationSolver = rlxslv;
  }


/** @} ---------------------------------------------------------------------*/
/*----------------- METHODS FOR ACCESSING THE DATA OF THE Block ------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the data of the Block
 *
 *
 *  @{ */


/** @} ---------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name 
 *  @{ */

 /// (try to) solve the Block

 int compute( bool changedvars = true ) override;

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessing the found solutions (if any) and solution information-
 *  @{ */

 void get_var_solution( Configuration * solc = nullptr ) override {
  f_incumbent->write( f_Block );  
 }

/*--------------------------------------------------------------------------*/
/// return the value of the (current) solution
/** Return the the value of the current solution. */

OFValue get_var_value() override { return( f_obj ); }

/** @} ---------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE Solver ----------------*/
/*--------------------------------------------------------------------------*/

 
/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of the BranchAndXSolver
 *
 *  @{ */
 
/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED FIELDS  -----------------------------*/
/*--------------------------------------------------------------------------*/
 
 Solver * f_RelaxationSolver;          // pointer to the relaxation solver
 RelaxationSolver * RlxSlv;
 Solution * f_incumbent;               // incumbent solution 

 double f_lb;                          // global lower bound
 double f_ub;                          // global upper bound
 double f_obj;                         // value of the objective function 

 int f_sense;                          // type of the optimization (min,max)
 
 // algorithmic parameters - - - - - - - - - - - - - - - - - - - - - - - - - 


/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// process all the pending modifications 

 void process_outstanding_Modification() {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// update current subproblem and its path and reverse path from the root

 void update_subproblem( BaXnode * old_node , BaXnode * node );
 
/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS  -----------------------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class BranchAndXSolver )

/** @} end( group( BranchAndXSolver_CLASSES ) ) ----------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* BranchAndXSolver.h included */

/*--------------------------------------------------------------------------*/
/*----------------------- End File BranchAndXSolver.h ----------------------*/
/*--------------------------------------------------------------------------*/
