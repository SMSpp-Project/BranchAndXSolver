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

class Node {
 
 public:
 
 using Index = Block::Index;

  // - - - - - - - - - - Constructor and Destructor - - - - - - - - - - - - -
 
  Node( Block * block = nullptr ): 
        f_block( block ) , 
        f_bound( + Inf< double >() ) {}
  
  ~Node() { delete f_block; };
 
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  void set_bound( double bound ) { f_bound = bound; }
  
  double bound() { return f_bound; }
  Block * block() { return f_block; }

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 protected:
  Block * f_block;
  double f_bound;

}; // end( class Node )

/*--------------------------------------------------------------------------*/

class Queue {
 public:
    Queue( Node * root ) { push( root ); }
    ~Queue() = default;

    bool empty() { return Q.empty(); }
    void push( Node * N ) { Q.push_back( N ); }
    Node * pop() { Node * N = Q.back(); Q.pop_back(); return N; }

    std::vector< Node * > Q;

}; // end( class Queue )

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
                      f_lb( - Inf< double >() ) ,
                      f_ub( + Inf< double >() ) ,
                      f_obj( + Inf< double >() ) ,
                      f_n_nodes( 0 ) ,  
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

 /*--------------------------------------------------------------------------*/ 

 Index get_n_nodes() { return f_n_nodes; }

/*--------------------------------------------------------------------------*/

 bool Max() { return( f_sense == Objective::eMax ); }

/*--------------------------------------------------------------------------*/

 bool Min() { return( f_sense == Objective::eMin ); }

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

 Index f_n_nodes;                      // number of nodes 
 
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
