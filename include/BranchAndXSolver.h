/*--------------------------------------------------------------------------*/
/*---------------------- File BranchAndXSolver.h -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class BranchAndXSolver, which implements
 * the Solver concept [see Solver.h] for a Relaxation-Agnostic Branch-and-X
 * (RABaX) Solver: the problem encoded by the Block is solved with an
 * enumerative algorithm built on top of the ChangeSolver / RelaxationSolver
 * concepts [see ChangeSolver.h]. The nodes of the enumeration tree are
 * identified by the Change that generates them from their parent, the
 * attached RelaxationSolver(s) provide dual bounds, true (primal) solutions
 * and the branching Changes, and the attached heuristic ChangeSolver(s), if
 * any, provide further primal bounds. The tree can be explored depth-first,
 * breadth-first or best-first (see intSolveMethod); nodes are pruned by the
 * standard inherited tolerances (dblRelAcc / dblAbsAcc), and intMaxNodes /
 * the inherited dblMaxTime bound the number of explored nodes and the total
 * time.
 *
 * The solver is generic ("X" is for Bound / Cut / Price): nothing in here
 * depends on the specific :Block being solved, all the problem-specific
 * knowledge lives in the attached :ChangeSolver, which provides the dual
 * bounds and the branching Changes (any tightening it does --- preprocessing,
 * reduced-cost fixing, cuts --- it does on its own terms, reading the search-
 * global data through a GlobalInformation [see ChangeSolver.h]). The tree is
 * explored serially or in parallel, both at the tree level (see
 * ParallelDFSSolve()) and across the solvers of a single node (see
 * intThreadForDifferentSolvers), and can be retained across re-solves for
 * reoptimization (see intReoptimize).
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
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __BranchAndXSolver
 #define __BranchAndXSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <atomic>
#include <chrono>
#include <climits>
#include <list>
#include <mutex>
#include <queue>

#include "BlockSolverConfig.h"

#include "ChangeSolver.h"

#include "Solution.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

 class Node;            // forward declaration of Node
 class DFSNode;         // forward declaration of DFSNode
 class ExploringNode;   // forward declaration of ExploringNode
 class OpenList;        // forward declaration of OpenList

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS BranchAndXSolver ------------------------*/
/*--------------------------------------------------------------------------*/
/// a generic Branch-and-Bound Solver driven by ChangeSolver
/** The BranchAndXSolver explores the enumeration tree whose nodes are
 * identified by the Change leading to them from their parent; moving the
 * attached :ChangeSolver between two nodes amounts to applying the Changes
 * (or the undo Changes) found along the joining path. RelaxationSolver(s)
 * provide the dual bound, possibly true solutions, and the branching
 * Changes; plain ChangeSolver(s) act as primal heuristics. */

class BranchAndXSolver : public Solver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 @{ */

 using Index = Block::Index;

 /// public enum for the int algorithmic parameters

 enum int_par_type_BXS {
  intSolveMethod = intLastAlgPar ,  ///< how to explore the tree
  intThreadForDifferentSolvers ,    ///< max threads for solvers at each node
  intReoptimize ,                   /**< retain the tree to reoptimize:
   * with it the (BestFS, serial) exploration keeps the tree and its fenced
   * frontier alive across compute() calls, and a re-solve under class 1-3
   * changes [see RelaxationSolver::classify()] re-seeds the search from
   * the re-evaluated frontier instead of re-deriving the whole tree. Note
   * that in a binary enumeration the fenced frontier is about as large as
   * the interior, so this pays off when evaluating a node is expensive
   * (say, LP-like relaxations) and/or few nodes re-open, while for very
   * cheap relaxations re-evaluating the frontier may cost as much as
   * solving from scratch; it is off by default, also because retaining
   * the tree can use a lot of memory. */
  intMaxNodes ,                     ///< max number of nodes to explore
  intLastBXSPar                     ///< first allowed new int parameter
  };

 /// public enum for the string algorithmic parameters

 enum str_par_type_BXS {
  /// name of the BlockSolverConfig file providing the inner Solver
  strNameOfBlockSolverConfigurationFile = strLastAlgPar ,
  strLastBXSPar                    ///< first allowed new string parameter
  };

 /// public enum for the possible tree exploration strategies

 enum SolveMethod {
  DFS = 0 ,                         ///< depth-first search
  BFS = 1 ,                         ///< breadth-first search
  BestFS = 2 ,                      ///< best-first search
  BestFSDive = 3                    ///< best-first search with depth-first dives
  };

/** @} ---------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor: the inner Solver are provided later via the configuration

 BranchAndXSolver() : Solver() , f_RelaxationSolvers() ,
                      f_HeuristicSolvers() , f_state( kUnEval ) ,
                      bestBound( 0 ) , bestSolution( nullptr ) ,
                      solveType( BestFS ) , maxThreadForSolvers( 1 ) ,
                      maxThread( 0 ) , reoptimize( 0 ) ,
                      maxNodes( INT_MAX ) , f_treeRoot( nullptr ) ,
                      nodeBudget( INT_MAX ) , timeBudget( Inf< double >() ) ,
                      relTol( 0 ) , absTol( 0 ) ,
                      configurationRS( nullptr ) , changes( 4 ) ,
                      subOptimalNodes() , infeasibleNodes() ,
                      integerNodes() {
  // the relaxations read the incumbent (for reduced-cost fixing and the like)
  // through the global information, kept bound to the live best-bound cell
  f_globalInfo.bind_incumbent( &bestBound );
  }

/*--------------------------------------------------------------------------*/
 /// constructor taking the inner Solver and the Block directly

 BranchAndXSolver( RelaxationSolver * Rsolver , ChangeSolver * Hsolver ,
                   Block * B ) : BranchAndXSolver() {
  if( Rsolver )
   f_RelaxationSolvers.push_back( Rsolver );
  if( Hsolver )
   f_HeuristicSolvers.push_back( Hsolver );
  Solver::set_Block( B );
  for( auto s : f_RelaxationSolvers )
   s->set_Block( B );
  for( auto s : f_HeuristicSolvers )
   s->set_Block( B );
  }

/*--------------------------------------------------------------------------*/
 /// destructor

 ~BranchAndXSolver() override {
  discardRetainedTree();
  delete bestSolution;
  delete configurationRS;
  for( auto slvr : v_created ) {
   slvr->set_Block( nullptr );
   delete slvr;
   }
  }

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// set the (pointer to the) Block that the Solver has to solve

 void set_Block( Block * block ) override {
  Solver::set_Block( block );
  if( f_RelaxationSolvers.empty() && f_HeuristicSolvers.empty() &&
      configurationRS )
   applyConfigurationToSolvers();
  }

/*--------------------------------------------------------------------------*/
 /// set the int parameters of BranchAndXSolver
 /** Set the int parameters specific of BranchAndXSolver:
  *
  * - intSolveMethod [BestFS]: how to explore the tree (see SolveMethod);
  *
  * - intThreadForDifferentSolvers [1]: maximum number of threads used to
  *   evaluate the (multiple) inner Solver at each node in parallel [see
  *   computeRelaxationsParallel() / computeHeuristicParallel()]; 1 = serial.
  *   It only pays when several (expensive) relaxations / heuristics are
  *   attached: their compute() is run concurrently and the results reduced
  *   exactly as in the serial path. Active in all the (serial-tree)
  *   explorations: depth-, breadth- and best-first.
  *
  * - intMaxNodes [no limit]: the budget of tree nodes a solve may explore;
  *   <= 0 means no limit. It is the enumerative counterpart of the inherited
  *   intMaxIter (a node is not an iteration: a node may take several inner
  *   iterations), which is left to the inner Solvers.
  *
  * The inherited intMaxThread [0] sets the number of workers of the
  * parallel tree exploration [see ParallelDFSSolve(); <= 1 = serial; the
  * base Solver classes do not store it, so it is stored here]. The inherited
  * dblMaxTime / dblRelAcc / dblAbsAcc respectively bound the total solve time
  * and set the optimality tolerances of the pruning. */

 void set_par( idx_type par , int value ) override {
  switch( par ) {
   case( intSolveMethod ):
    solveType = static_cast< SolveMethod >( value );
    break;
   case( intMaxThread ):
    maxThread = std::max( 0 , value );
    break;
   case( intThreadForDifferentSolvers ):
    if( value <= 0 )
     throw( std::invalid_argument( "BranchAndXSolver::set_par: "
            "intThreadForDifferentSolvers must be positive" ) );
    maxThreadForSolvers = value;
    break;
   case( intReoptimize ):
    reoptimize = value;
    break;
   case( intMaxNodes ):
    maxNodes = ( value <= 0 ) ? INT_MAX : value;
    break;
   default:
    Solver::set_par( par , value );
   }
  }

/*--------------------------------------------------------------------------*/
 /// set the string parameters of BranchAndXSolver
 /** Set the string parameters specific of BranchAndXSolver:
  *
  * - strNameOfBlockSolverConfigurationFile [""]: name of the file out of
  *   which the BlockSolverConfig providing the inner Solver (the
  *   RelaxationSolver(s) and heuristic ChangeSolver(s)) is deserialized;
  *   it is applied to the Block as soon as both are available. */

 void set_par( idx_type par , const std::string & value ) override {
  if( par != strNameOfBlockSolverConfigurationFile ) {
   Solver::set_par( par , value );
   return;
   }
  nameRS = value;
  configurationRS = dynamic_cast< BlockSolverConfig * >(
                                      Configuration::deserialize( nameRS ) );
  if( ! configurationRS )
   throw( std::invalid_argument( "BranchAndXSolver::set_par: unable to "
          "create a BlockSolverConfig out of " + nameRS ) );
  if( f_Block )
   applyConfigurationToSolvers();
  }

/*--------------------------------------------------------------------------*/
 /// receive a Modification and forward it to the (private) inner Solver

 void add_Modification( sp_Mod & mod ) override {
  Solver::add_Modification( mod );
  for( auto slvr : v_created )
   slvr->add_Modification( mod );
  }

/** @} ---------------------------------------------------------------------*/
/*--------------------- METHODS FOR SOLVING THE MODEL ----------------------*/
/*--------------------------------------------------------------------------*/

 /// solve the problem encoded in the Block with Branch-and-Bound

 int compute( bool changedvars = true ) override;

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

 /// write the best solution found into the Variable of the Block

 void get_var_solution( Configuration * solc = nullptr ) override {
  if( bestSolution )
   bestSolution->write( f_Block );
  }

/*--------------------------------------------------------------------------*/
 /// return the value of the best solution found

 OFValue get_var_value( void ) override { return( bestBound ); }

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR HANDLING THE PARAMETERS ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the parameters of the BranchAndXSolver
 *  @{ */

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( idx_type( intLastBXSPar ) );
  }

 [[nodiscard]] idx_type get_num_str_par( void ) const override {
  return( idx_type( strLastBXSPar ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override {
  switch( par ) {
   case( intSolveMethod ):               return( int( BestFS ) );
   case( intThreadForDifferentSolvers ): return( 1 );
   case( intReoptimize ):                return( 0 );
   case( intMaxNodes ):                  return( INT_MAX );
   }
  return( Solver::get_dflt_int_par( par ) );
  }

 [[nodiscard]] const std::string & get_dflt_str_par( idx_type par )
  const override {
  if( par == strNameOfBlockSolverConfigurationFile ) {
   static const std::string dfltRS;
   return( dfltRS );
   }
  return( Solver::get_dflt_str_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] int get_int_par( idx_type par ) const override {
  switch( par ) {
   case( intSolveMethod ):               return( int( solveType ) );
   case( intMaxThread ):                 return( maxThread );
   case( intThreadForDifferentSolvers ): return( maxThreadForSolvers );
   case( intReoptimize ):                return( reoptimize );
   case( intMaxNodes ):                  return( maxNodes );
   }
  return( Solver::get_int_par( par ) );
  }

 [[nodiscard]] const std::string & get_str_par( idx_type par )
  const override {
  if( par == strNameOfBlockSolverConfigurationFile )
   return( nameRS );
  return( Solver::get_str_par( par ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override {
  if( name == "intSolveMethod" )
   return( intSolveMethod );
  if( name == "intThreadForDifferentSolvers" )
   return( intThreadForDifferentSolvers );
  if( name == "intReoptimize" )
   return( intReoptimize );
  if( name == "intMaxNodes" )
   return( intMaxNodes );
  return( Solver::int_par_str2idx( name ) );
  }

 [[nodiscard]] idx_type str_par_str2idx( const std::string & name )
  const override {
  if( name == "strNameOfBlockSolverConfigurationFile" )
   return( strNameOfBlockSolverConfigurationFile );
  return( Solver::str_par_str2idx( name ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override {
  static const std::string pars[] = { "intSolveMethod" ,
                                      "intThreadForDifferentSolvers" ,
                                      "intReoptimize" , "intMaxNodes" };
  if( ( idx >= intSolveMethod ) && ( idx < intLastBXSPar ) )
   return( pars[ idx - intSolveMethod ] );
  return( Solver::int_par_idx2str( idx ) );
  }

 [[nodiscard]] const std::string & str_par_idx2str( idx_type idx )
  const override {
  if( idx == strNameOfBlockSolverConfigurationFile ) {
   static const std::string par = "strNameOfBlockSolverConfigurationFile";
   return( par );
   }
  return( Solver::str_par_idx2str( idx ) );
  }

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED METHODS ----------------------------*/
/*--------------------------------------------------------------------------*/

 /// setter for the tree exploration strategy

 void setSolveType( SolveMethod st ) { solveType = st; }

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS ----------------------------*/
/*--------------------------------------------------------------------------*/

 /// pointer(s) to the Solver used to solve the relaxations
 std::vector< RelaxationSolver * > f_RelaxationSolvers;

 /// pointer(s) to the Solver used to find feasible solutions
 std::vector< ChangeSolver * > f_HeuristicSolvers;

 int f_state;             ///< the (current) state of the compute() process

 double bestBound;        ///< the best primal (hence feasible) bound found

 Solution * bestSolution; ///< the best solution found, valued bestBound

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

 /// create the inner Solver out of configurationRS
 /** Creates and ComputeConfig-ures the Solver described by the
  * BlockSolverConfig and attaches them to the Block, WITHOUT registering
  * them: the inner Solver are private to the BranchAndXSolver, which
  * owns them, drives them, and forwards them the Modification it receives
  * [see add_Modification()]; this way they never appear in the Block's
  * registered-Solver list (whose order and content belong to the user). */

 /// create the per-worker private Solver sets for the parallel exploration
 /** Creates \p K additional sets of inner Solver out of configurationRS,
  * one per worker of the parallel exploration (see v_workerSolvers): each
  * worker drives its own clones, so that the (cheap, internal) application
  * of the branching Changes never needs synchronization; the incumbent and
  * the pool of open subtrees are the only shared state. The clones are
  * attached to the Block (set_Block, no registration) and owned. */

 void createWorkerSolvers( Index K ) {
  if( ! ( configurationRS && f_Block ) )
   throw( std::logic_error( "BranchAndXSolver::createWorkerSolvers: "
          "configurationRS or f_Block not set" ) );
  v_workerSolvers.resize( K );
  for( Index w = 0 ; w < K ; ++w ) {
   auto & ws = v_workerSolvers[ w ];
   if( ! ws.relaxation.empty() )    // already created (re-solve)
    continue;
   for( Index i = 0 ; i < configurationRS->num_ComputeConfig() ; ++i ) {
    auto slvr = Solver::new_Solver( configurationRS->get_SolverName( i ) );
    if( auto cfg = configurationRS->get_SolverConfig( i ) )
     slvr->set_ComputeConfig( cfg );
    slvr->set_Block( f_Block );
    v_created.push_back( slvr );    // owned like the serial ones
    if( auto rs = dynamic_cast< RelaxationSolver * >( slvr ) ) {
     rs->set_global_information( &f_globalInfo );
     ws.relaxation.push_back( rs );
     }
    else if( auto hs = dynamic_cast< ChangeSolver * >( slvr ) )
     ws.heuristic.push_back( hs );
    else
     throw( std::invalid_argument( "BranchAndXSolver::"
            "createWorkerSolvers: the BlockSolverConfig must only "
            "describe :ChangeSolver" ) );
    }
   }
  }

/*--------------------------------------------------------------------------*/

 void applyConfigurationToSolvers( void ) {
  if( ! ( configurationRS && f_Block ) )
   throw( std::logic_error( "BranchAndXSolver::"
          "applyConfigurationToSolvers: configurationRS or f_Block not set"
          ) );
  for( Index i = 0 ; i < configurationRS->num_ComputeConfig() ; ++i ) {
   auto slvr = Solver::new_Solver( configurationRS->get_SolverName( i ) );
   if( auto cfg = configurationRS->get_SolverConfig( i ) )
    slvr->set_ComputeConfig( cfg );
   slvr->set_Block( f_Block );
   v_created.push_back( slvr );
   if( auto rs = dynamic_cast< RelaxationSolver * >( slvr ) ) {
    rs->set_global_information( &f_globalInfo );
    f_RelaxationSolvers.push_back( rs );
    }
   else if( auto hs = dynamic_cast< ChangeSolver * >( slvr ) )
    f_HeuristicSolvers.push_back( hs );
   else
    throw( std::invalid_argument( "BranchAndXSolver::"
           "applyConfigurationToSolvers: the BlockSolverConfig must only "
           "describe :ChangeSolver" ) );
   }
  }

/*--------------------------------------------------------------------------*/
 /// parallel depth-first exploration with \p K workers
 /** Parallel depth-first exploration, used when the inherited intMaxThread
  * parameter is > 1: a serial, ordered (FIFO, i.e., discovery order)
  * ramp-up expands the tree BestFS-style until enough open subtrees exist,
  * then \p K workers, each driving its own private set of inner Solver
  * [see createWorkerSolvers()], repeatedly claim the OLDEST open subtree
  * (preserving the sequential search order, which is what keeps parallel
  * performance replicable) and explore it depth-first, sharing only the
  * incumbent (under mutex) and the pool of open subtrees.
  *
  * The distribution of work is dynamic: when the pool runs low a worker
  * DONATES one of the (eagerly evaluated, hence claimable) top-level
  * branches of its current subtree instead of exploring it, so that
  * deep-and-narrow trees keep all the workers busy; donated children hang
  * off the claimed node, which always outlives them, keeping their path
  * replayable and the final cleanup single-owner. Plain std::thread is
  * used on purpose: at this granularity (one whole subtree per claim) a
  * task framework would add nothing. */

 int ParallelDFSSolve( std::mutex & globalMutex , int K );

/*--------------------------------------------------------------------------*/
 /// delete the tree retained for reoptimization (if any) and its frontier

 void discardRetainedTree( void );

/*--------------------------------------------------------------------------*/
 /// the depth-first subtree exploration of one worker
 /** The recursive core of one worker of ParallelDFSSolve(): a depth-first
  * exploration on the worker's own Solver set, with the shared
  * incumbent updated under \p incumbentMutex, a shared atomic node budget
  * and a wall-clock deadline. */

 int workerDFS( Node * currentNode , std::list< ChangeSolver * > & solvers ,
                std::vector< RelaxationSolver * > & relaxation ,
                std::vector< ChangeSolver * > & heuristic ,
                bool minimizing , std::mutex & incumbentMutex ,
                std::atomic< int > & nodeBdg ,
                std::chrono::steady_clock::time_point deadline ,
                int & nameCounter );

/*--------------------------------------------------------------------------*/
 /// solve the tree breadth-first
 /** @param globalMutex mutex protecting the shared state of this class
  *  @return the sol_type [see Solver.h] of the computation */

 int BFSSolve( std::mutex & globalMutex );

/*--------------------------------------------------------------------------*/
 /// solve the tree depth-first
 /** Depth-first exploration: the same explore() loop as the breadth- and
  * best-first ones, only with a LIFO stack as the open set, over the common
  * ExploringNode tree.
  *  @param globalMutex mutex protecting the shared state of this class
  *  @return the sol_type [see Solver.h] of the computation */

 int DFSSolve( std::mutex & globalMutex );

/*--------------------------------------------------------------------------*/
 /// solve the tree best-first
 /** @param globalMutex mutex protecting the shared state of this class
  *  @return the sol_type [see Solver.h] of the computation */

 int BestFirstSolve( std::mutex & globalMutex );

/*--------------------------------------------------------------------------*/
 /// the single serial tree exploration, driven by the open-list discipline
 /** The one exploration loop shared by every (serial) strategy: it
  * repeatedly takes the next node from \p open [see OpenList], moves the
  * :ChangeSolver to it, evaluates and branches it, and pushes the surviving
  * children back into \p open; the discipline of \p open (stack / queue /
  * priority) is what makes it a depth- / breadth- / best-first search.
  * \p currentNode is where the :ChangeSolver currently sits (the root, or
  * the last re-seeded node in a reoptimization), \p retain keeps the tree
  * and its fenced frontier alive for a later reoptimization [see
  * intReoptimize]. */

 int explore( OpenList & open , std::mutex & globalMutex ,
              std::list< ChangeSolver * > * solvers , bool minimizing ,
              int & counter , ExploringNode * rootNode ,
              ExploringNode * currentNode , bool retain ,
              std::chrono::high_resolution_clock::time_point start );

/*--------------------------------------------------------------------------*/
 /// tear down a retained / in-progress tree and its open set and frontier

 void discardTree( ExploringNode * root , OpenList & open ,
                   std::list< ChangeSolver * > * solvers );

/*--------------------------------------------------------------------------*/
 /// process the outstanding Modification
 /** Processes the queued Modification and classifies the kind of change
  * they amount to (see the changes field), in view of reoptimizing the
  * tree rather than rebuilding it; the reoptimization itself is NOT
  * implemented yet, any change currently triggers a solve from scratch. */

 void process_outstanding_Modification( void );

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 SolveMethod solveType;       ///< chosen method to explore the tree

 int maxThreadForSolvers;     ///< max threads for solvers at each node

 int maxThread;               ///< workers of the parallel tree exploration

 int reoptimize;              ///< retain the tree to reoptimize (see
                              ///< intReoptimize / BestFirstSolve())

 int maxNodes;                ///< node budget of a solve (see intMaxNodes)

 /// the search-global information shared with the relaxations (incumbent,
 /// global cuts/columns); its incumbent is bound to the live bestBound
 GlobalInformation f_globalInfo;

 /// the root of the tree retained for reoptimization, nullptr if none
 ExploringNode * f_treeRoot;

 /// residual node budget of the current solve (from intMaxNodes), consumed
 /// by the exploration
 int nodeBudget;

 /// residual time budget of the current solve (from the inherited
 /// dblMaxTime), consumed by the exploration
 double timeBudget;

 double relTol;   ///< the inherited dblRelAcc, cached at compute() start

 double absTol;   ///< the inherited dblAbsAcc, cached at compute() start

 std::string nameRS;          ///< name of the BlockSolverConfig file

 BlockSolverConfig * configurationRS;  ///< the inner BlockSolverConfig

 /// the (private) inner Solver created out of configurationRS, owned
 std::vector< Solver * > v_created;

 /// the private Solver set of one worker of the parallel exploration
 struct WorkerSolvers {
  std::vector< RelaxationSolver * > relaxation;
  std::vector< ChangeSolver * > heuristic;
  };

 /// the per-worker Solver sets [see createWorkerSolvers()]
 std::vector< WorkerSolvers > v_workerSolvers;

 /// classification of the outstanding changes, for future reoptimization:
 /// 0 = none, 1 = objective only, 2 = r.h.s. only, 3 = both objective and
 /// r.h.s., anything else = general change (solve from scratch)
 char changes;

 /// nodes pruned as sub-optimal, to be revisited on an objective change
 std::list< ExploringNode * > subOptimalNodes;

 /// nodes pruned as infeasible, to be revisited on a r.h.s. change
 std::list< ExploringNode * > infeasibleNodes;

 /// nodes pruned as integer-feasible, to be revisited on reoptimization
 std::list< ExploringNode * > integerNodes;

 SMSpp_insert_in_factory_h;   // insert BranchAndXSolver in the factory

/*--------------------------------------------------------------------------*/

 };  // end( class( BranchAndXSolver ) )

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASS Node --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a node of the Branch-and-Bound tree
/** A Node of the enumeration tree is identified by the Change that leads to
 * it from its parent; it also holds the (undo) Change to return to the
 * parent, the Changes generating its children (see obtainBranchList()) and
 * its dual bound. The Node owns all these Change. */

class Node {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
 /// constructor
 /** @param change the Change to apply to reach this node from its parent
  *  @param nodeName name of the node (debugging purposes only) */

 Node( Change * change , int nodeName = 0 )
  : dual_bound( - Inf< double >() ) , f_change( change ) ,
    toFather( nullptr ) , branches() , f_infeasible( false ) ,
    name( nodeName ) {}

/*--------------------------------------------------------------------------*/
 /// destructor: deletes the owned Changes

 virtual ~Node() {
  delete f_change;
  delete toFather;
  for( auto ch : branches )
   delete ch;
  }

/*--------------------------------------------------------------------------*/

 /// return the value of the dual bound
 double get_dual_bound( void ) const { return( dual_bound ); }

 /// set the value of the dual bound
 void set_dual_bound( double db ) { dual_bound = db; }

 /// return the Change that leads to this node
 Change * get_f_change( void ) const { return( f_change ); }

 /// return the Change to undo to return to the parent node
 Change * get_toFather( void ) const { return( toFather ); }

 /// set the Change to undo to return to the parent node
 void set_toFather( Change * tf ) { toFather = tf; }

 /// return the name of the node (debugging purposes only)
 int get_name( void ) const { return( name ); }

 /// initialize the dual bound to the appropriate infinity
 void initializeBound( bool minimizing ) {
  dual_bound = minimizing ? - Inf< double >() : Inf< double >();
  f_infeasible = false;
  }

 /// whether the node was found infeasible (its relaxation has no solution)
 bool is_infeasible( void ) const { return( f_infeasible ); }

 /// record that the node is infeasible
 void set_infeasible( bool i ) { f_infeasible = i; }

 /// have the given RelaxationSolver produce the branching Changes
 void obtainBranchList( RelaxationSolver * solver ) {
  branches = solver->branch();
  }

 /// return the Changes generating the children of this node
 std::vector< Change * > & getBranches( void ) { return( branches ); }

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

 double dual_bound;     ///< value of the dual bound at this node

 Change * f_change;     ///< the Change to reach this node from the parent

 Change * toFather;     ///< the (undo) Change to return to the parent

 /// the Changes to reach each child from this node
 std::vector< Change * > branches;

 bool f_infeasible;     ///< whether the node's relaxation has no solution

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 int name;              ///< name of the node, debugging purposes only

/*--------------------------------------------------------------------------*/

 };  // end( class( Node ) )

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASS DFSNode ------------------------------*/
/*--------------------------------------------------------------------------*/
/// the Node used in the depth-first exploration of the tree
/** In the recursive depth-first exploration no explicit tree structure is
 * needed (the recursion stack is the path to the root), so a DFSNode is
 * just a plain Node. */

class DFSNode : public Node {

 public:

 /// constructor, see Node::Node()

 DFSNode( Change * change , int nodeName = 0 ) : Node( change , nodeName ) {}

 ~DFSNode() override = default;   ///< destructor

 };  // end( class( DFSNode ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS ExploringNode ---------------------------*/
/*--------------------------------------------------------------------------*/
/// the Node used in the breadth-first / best-first explorations
/** An ExploringNode also knows its parent and children, so that the
 * :ChangeSolver can be moved between two arbitrary nodes of the tree (see
 * moveBetweenNodes()) by climbing up and down along the joining path. */

class ExploringNode : public Node {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
 /// constructor
 /** @param change the Change to apply to reach this node from its parent
  *  @param parent pointer to the parent node
  *  @param depth depth of the node in the tree
  *  @param nodeName name of the node (debugging purposes only) */

 ExploringNode( Change * change , ExploringNode * parent , int depth ,
                int nodeName = 0 )
  : Node( change , nodeName ) , f_children() , level( depth ) ,
    parent( parent ) {}

 ~ExploringNode() override = default;   ///< destructor

/*--------------------------------------------------------------------------*/

 /// getter for the level of the node
 int get_level( void ) const { return( level ); }

 /// getter for the parent node
 ExploringNode * get_parent( void ) const { return( parent ); }

 /// setter for the parent node
 void set_parent( ExploringNode * p ) { parent = p; }

 /// getter for the list of children nodes
 std::list< ExploringNode * > & get_children( void ) {
  return( f_children );
  }

/*--------------------------------------------------------------------------*/
 /// move the given :ChangeSolver between two nodes of the tree
 /** Moves the given :ChangeSolver from \p sourceNode to \p destNode by
  * applying, along the path joining them, the undo Changes (climbing up
  * from the source to the common ancestor) and then, top-down, the node
  * Changes (descending from the common ancestor to the destination). */

 static void moveBetweenNodes( ExploringNode * sourceNode ,
                               ExploringNode * destNode ,
                               std::list< ChangeSolver * > * solvers ) {
  std::list< Change * > changesToApply;
  while( sourceNode->level != destNode->level ) {
   if( sourceNode->level > destNode->level ) {
    for( const auto s : *solvers )
     s->apply( sourceNode->toFather , false );
    sourceNode = sourceNode->parent;
    }
   else {
    changesToApply.push_front( destNode->f_change );
    destNode = destNode->parent;
    }
   }
  while( sourceNode != destNode ) {
   for( const auto s : *solvers )
    s->apply( sourceNode->toFather , false );
   sourceNode = sourceNode->parent;
   changesToApply.push_front( destNode->f_change );
   destNode = destNode->parent;
   }
  for( auto change : changesToApply )
   for( const auto s : *solvers )
    s->apply( change , false );
  }

/*--------------------------------------------------------------------------*/
 /// prune the node, recursively pruning any parent left childless
 /** Prunes \p nodeToPrune, moving the :ChangeSolver to its parent; if the
  * parent is left with no children it is recursively pruned as well.
  *  @return the node the :ChangeSolver is left at, i.e., the parent of the
  *          last pruned node */

 static ExploringNode * prune( ExploringNode * nodeToPrune ,
                               std::list< ChangeSolver * > * solvers ) {
  ExploringNode * parentNode = nodeToPrune->get_parent();
  if( nodeToPrune->get_toFather() )
   for( const auto s : *solvers )
    s->apply( nodeToPrune->get_toFather() , false );
  if( ! parentNode->get_children().empty() )
   parentNode->get_children().remove( nodeToPrune );
  auto & branches = parentNode->getBranches();
  branches.erase( std::remove( branches.begin() , branches.end() ,
                               nodeToPrune->get_f_change() ) ,
                  branches.end() );
  delete nodeToPrune;
  if( parentNode->get_children().empty() && parentNode->get_parent() )
   return( prune( parentNode , solvers ) );
  return( parentNode );
  }

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

 std::list< ExploringNode * > f_children;  ///< list of the children nodes

 int level;                                ///< level of the node in the tree

 ExploringNode * parent;                   ///< pointer to the parent node

/*--------------------------------------------------------------------------*/

 };  // end( class( ExploringNode ) )

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASS OpenList ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the set of open nodes of an exploration, abstracted over its discipline
/** The frontier of open nodes that a tree exploration keeps and picks the
 * next node from. Abstracting the discipline (which node comes next) behind
 * this interface lets a single exploration loop [see
 * BranchAndXSolver::explore()] serve every strategy: a LIFO stack gives
 * depth-first, a FIFO queue gives breadth-first, a priority queue ordered by
 * dual bound gives best-first. A new strategy is just a new OpenList. */

class OpenList {

 public:

 virtual ~OpenList() = default;

 /// whether there are no open nodes left
 [[nodiscard]] virtual bool empty( void ) const = 0;

 /// add a node to the open set
 virtual void push( ExploringNode * node ) = 0;

 /// remove and return the next node to explore
 virtual ExploringNode * pop( void ) = 0;

 /// whether the discipline is last-in first-out (a stack)
 /** Tells the exploration how to order a node's freshly generated children
  * in the open set so that the most promising one [the first returned by
  * RelaxationSolver::branch()] is explored first: a LIFO stack must receive
  * them in reverse, the other disciplines in branching order. */

 [[nodiscard]] virtual bool isLIFO( void ) const { return( false ); }

 };  // end( class( OpenList ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* BranchAndXSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------------- End File BranchAndXSolver.h ----------------------*/
/*--------------------------------------------------------------------------*/
