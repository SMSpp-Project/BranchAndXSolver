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

#include <ctime>

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

 BranchAndXSolver() : f_Solver( nullptr ) ,
                      f_RelaxationSolver( nullptr ) ,
                      f_incumbent( nullptr ) , 
                      f_lb( - Inf< double >() ) ,
                      f_ub( + Inf< double >() ) ,
                      f_n_nodes( 0 ) ,  
                      f_sense( Objective::eMin ) {

 // Initiliaze the value of the objective function
 f_obj =  f_sense == Objective::eMin ? + Inf< double >() : - Inf< double >();

 // Algorithmic parameters: default values
 MaxIter = 100000;
 MaxTime = 3600;
 RelAcc = 1e-04;
 AbsAcc = 1e-04;
 }

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

 void set_RelaxationSolver( Solver * slv ) {
  
  f_Solver = slv;  
  f_RelaxationSolver = dynamic_cast< RelaxationSolver * >( slv );
  if( ! f_RelaxationSolver )
   throw( std::invalid_argument( "Not able to cast to RelaxationSolver *" ) );

  }

/*--------------------------------------------------------------------------*/ 

 Index get_n_nodes() { return f_n_nodes; }

/*--------------------------------------------------------------------------*/

 int get_objective_sense() { return f_sense; }

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

class Node {
 
 public:
 
 using Index = Block::Index;

  // - - - - - - - - - - Constructor and Destructor - - - - - - - - - - - - -
 
  Node( Block * block = nullptr ): 
        f_block( block ) , 
        f_ub( + Inf< double >() ) ,
        f_lb( - Inf< double >() ) {}
  
  ~Node() { delete f_block; };
 
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  bool is_feasible( double Acc ) { return ( std::abs( f_ub - f_lb ) < Acc ); }

  void set_ub( double ub ) { f_ub = ub; }
  void set_lb( double lb ) { f_lb = lb; }
  double get_ub() { return f_ub; }
  double get_lb() { return f_lb; }

  Block * block() { return f_block; }

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 protected:
  Block * f_block;
  double f_ub;
  double f_lb;

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
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED FIELDS  -----------------------------*/
/*--------------------------------------------------------------------------*/
 
 Solver * f_Solver;                     // pointer to the Solver
 RelaxationSolver * f_RelaxationSolver; // pointer to the Relaxation Solver
 Solution * f_incumbent;                // incumbent solution 

 Index f_n_nodes;                       // number of nodes 
 
 double f_lb;                           // global lower bound
 double f_ub;                           // global upper bound
 double f_obj;                          // value of the objective function 

 int f_sense;                           // type of the optimization (min,max)
 
 // algorithmic parameters - - - - - - - - - - - - - - - - - - - - - - - - - 

 double RelAcc;     ///< relative accuracy for declaring a solution optimal
 double AbsAcc;     ///< absolute accuracy for declaring a solution optimal
 double MaxTime;    ///< maximum time (in seconds) 
 Index MaxIter;     ///< maximum number of iterations

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
