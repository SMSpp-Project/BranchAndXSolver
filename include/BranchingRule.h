/*--------------------------------------------------------------------------*/
/*------------------------- File BranchingRule.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the components that a RelaxationSolver [see
 * ChangeSolver.h] can use to choose the variable to branch upon: the
 * description of a candidate and its features (BranchingFeatures), strong
 * branching with the product score (StrongBranching), and a scoring function
 * of the candidates read from a file (LearnedBranching); and, for the order
 * in which the children are explored, the estimates of the best value under
 * each of them (BranchEstimator), the pseudocosts that give them
 * (PseudoCosts), and the pseudocosts read from a file for the candidates
 * that have none yet (LearnedPseudoCosts). They live beside the
 * BranchAndXSolver and, but for BranchEstimator, know nothing of it: each
 * RelaxationSolver describes its own candidates, decides which component to
 * use (say, by a parameter of its own), and gives StrongBranching the way to
 * evaluate the children of a candidate, so that the same rule serves any
 * relaxation without the RelaxationSolver exposing anything of its
 * representation.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __BranchingRule
 #define __BranchingRule  /* self-identification: #endif at the end of file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <array>

#include <cstddef>

#include <functional>

#include <istream>

#include <ostream>

#include <string>

#include <unordered_map>

#include <vector>

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup BranchingRule_CLASSES Classes in BranchingRule.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*------------------------- STRUCT BranchCandidate -------------------------*/
/*--------------------------------------------------------------------------*/
/// a candidate for branching, as a RelaxationSolver describes it
/** A variable that the RelaxationSolver may branch upon at the current node,
 * described in terms that do not depend on the relaxation: the
 * RelaxationSolver has a point whose fractional components are the
 * candidates, a second point (which may be the same), and a partition of
 * the variables into blocks, each with its objective and with its
 * variables in a given order. The variable itself is only known to the
 * RelaxationSolver, by the index id. */

struct BranchCandidate {
 std::size_t id;      ///< the index of the variable, for the RelaxationSolver
 double value;        ///< its value in the point of the candidates
 double second;       ///< its value in the second point
 double cost;         ///< its cost in the objective of its block
 double position;     ///< its position among those of its block, in [0, 1]
 unsigned int block;  ///< the index of its block
 };

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS BranchingFeatures -------------------------*/
/*--------------------------------------------------------------------------*/
/// the fractionality of the candidates and their features
/** The most fractional candidate is the first one after
 * sort_by_fractionality(); the features of the candidates [see features()]
 * are what LearnedBranching scores. */

class BranchingFeatures
{
 public:

 /// the number of the features of a candidate [see features()]
 static constexpr std::size_t NFeatures = 8;

 /// the features of a candidate [see features()]
 using Features = std::array< double , NFeatures >;

 /// the fractionality of a value, i.e., its distance from the nearest integer
 [[nodiscard]] static double fractionality( double value );

 /// sort the candidates by nonincreasing fractionality
 /** The sort is stable, so that candidates with the same fractionality keep
  * the order in which the RelaxationSolver gave them. */

 static void sort_by_fractionality( std::vector< BranchCandidate > & cands );

 /// the features of the candidates
 /** In this order: (i) the fractionality of the candidate, (ii) its value
  * in the point of the candidates, (iii) its value in the second point,
  * (iv) the absolute difference of the two, (v) its cost over the largest
  * absolute cost of the candidates, (vi) its position in its block, (vii)
  * the fraction \p fixed of the variables that are fixed at the node and
  * (viii) its rank in \p cands over their number less one. They are
  * therefore relative to the candidates given, which are supposed to be
  * sorted by sort_by_fractionality(), and none of them depends on the
  * scale of the data. */

 [[nodiscard]] static std::vector< Features > features(
				const std::vector< BranchCandidate > & cands ,
				double fixed );

 };  // end( class( BranchingFeatures ) )

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS StrongBranching --------------------------*/
/*--------------------------------------------------------------------------*/
/// strong branching with the product score
/** For each candidate the bounds \f$ L^0_v \f$ and \f$ L^1_v \f$ of its two
 * children are computed by the function that the RelaxationSolver gives
 * [see Evaluator], and the score is the product
 * \f[
 *  s_v = \max \{ L^0_v - L , \delta \} \cdot \max \{ L^1_v - L , \delta \}
 * \f]
 * of the improvements over the bound \f$ L \f$ of the node (with the sign
 * of the sense of the problem), with \f$ \delta = 10^{-6} \max \{ | L | , 1
 * \} \f$, a child found empty counting as an infinite improvement. */

class StrongBranching
{
 public:

 /// what strong branching finds of a candidate: the bound of the "down"
 /// child, that of the "up" one and the score
 using Eval = std::array< double , 3 >;

 /// the bounds of the "down" and of the "up" child of a candidate
 /** Given by the RelaxationSolver, on the side of the relaxation (the lower
  * one for a minimization problem), +/- infinity for a child found empty;
  * the RelaxationSolver is expected to be back at the node afterwards. */

 using Evaluator =
  std::function< std::array< double , 2 >( const BranchCandidate & ) >;

 /// the bounds of the children and the score of each candidate
 /** \p node is the bound of the node and \p maximize tells the sense of the
  * problem; if \p log is not nullptr, one line per candidate is written to
  * it. */

 [[nodiscard]] static std::vector< Eval > evaluate(
				const std::vector< BranchCandidate > & cands ,
				const Evaluator & children , double node ,
				bool maximize , std::ostream * log = nullptr );

 /// the index of the candidate with the largest score
 /** Ties go to the first one, i.e., to the most fractional if the
  * candidates are sorted by BranchingFeatures::sort_by_fractionality(). */

 [[nodiscard]] static std::size_t best( const std::vector< Eval > & evals );

 /// append the data of strong branching to a file, in CSV
 /** One line per candidate: \p call (the progressive number of the
  * branching), \p node (the bound of the node), the fraction of the fixed
  * variables, the rank of the candidate, its block, its position, its
  * value in the point of the candidates and in the second one, its
  * fractionality, its cost, the bounds of the two children, the score, and
  * 1 for the candidate \p chosen; a header line is written when the file
  * is created, so that LearnedBranching can be fitted on these data.
  * Throws if the file cannot be opened. */

 static void write_log( const std::string & file , unsigned long call ,
			double node ,
			const std::vector< BranchCandidate > & cands ,
			const std::vector< BranchingFeatures::Features > & x ,
			const std::vector< Eval > & evals , std::size_t chosen );

 };  // end( class( StrongBranching ) )

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS LearnedBranching -------------------------*/
/*--------------------------------------------------------------------------*/
/// a scoring function of the candidates read from a file
/** A perceptron \f$ \sigma( \phi_v ) = w_2^\top \tanh( W_1 \phi_v + b_1 ) +
 * b_2 \f$ with one hidden layer, which gives a score to a candidate out of
 * its features \f$ \phi_v \f$ [see BranchingFeatures::features()], written
 * in a text file as the line "LBRModel <nin> <nhid>", with nin ==
 * BranchingFeatures::NFeatures, followed by the nhid rows of nin weights of
 * the hidden layer, its nhid biases, the nhid weights of the output and its
 * bias. */

class LearnedBranching
{
 public:

 /// read the model from the given file; throws if it is not one
 void load( const std::string & file );

 /// read the model from the given stream, named name in the errors; throws
 /// if it is not one, and leaves the stream after the model
 void load( std::istream & in , const std::string & name );

 /// whether a model has been read
 [[nodiscard]] bool empty( void ) const { return( f_w1.empty() ); }

 /// the score of a candidate with the given features
 [[nodiscard]] double score( const BranchingFeatures::Features & x ) const;

 /// the index of the candidate the model scores best
 /** Ties go to the first one. */

 [[nodiscard]] std::size_t best(
		     const std::vector< BranchingFeatures::Features > & x ) const;

 private:

 std::vector< std::vector< double > > f_w1;  ///< weights of the hidden layer
 std::vector< double > f_b1;                 ///< biases of the hidden layer
 std::vector< double > f_w2;                 ///< weights of the output
 double f_b2 = 0;                            ///< bias of the output

 };  // end( class( LearnedBranching ) )

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS BranchEstimator --------------------------*/
/*--------------------------------------------------------------------------*/
/// the estimates of the best value under the children of a branching
/** A RelaxationSolver [see ChangeSolver.h] that also derives from
 * BranchEstimator gives, right after each branch() and before any other of
 * its methods is called, an estimate of the value of the best solution of
 * the original problem in the region of each child, which the
 * BranchAndXSolver uses to order the open nodes [see intSolveMethod]. An
 * estimate is no bound: the nodes are still pruned by their bounds, so that
 * a wrong estimate only changes the order of the exploration. */

class BranchEstimator
{
 public:

 virtual ~BranchEstimator() = default;   ///< destructor: does nothing

 /// the estimates of the children of the last branch()
 /** One per Change returned by the last branch(), in the same order, in
  * the sense of the problem (an estimate of a minimum is expected to be no
  * less than the bound of the child); NaN for a child without one, and an
  * empty vector if there is none at all. */

 [[nodiscard]] virtual std::vector< double > branch_estimates( void ) = 0;

 };  // end( class( BranchEstimator ) )

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS PseudoCosts ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the pseudocosts of the candidates, and the estimates they give
/** The pseudocost \f$ \Psi^d_v \f$ of the candidate v in the direction d (0
 * "down", 1 "up") is the mean, over the branchings on v observed so far, of
 * the worsening of the bound from the node to its child in that direction
 * per unit of the distance moved by the value of v, which is \f$ x_v -
 * \lfloor x_v \rfloor \f$ down and \f$ \lceil x_v \rceil - x_v \f$ up. A
 * candidate never branched upon in a direction takes its prior, if one is
 * given, and otherwise the mean of the pseudocosts of all the candidates in
 * that direction (0 if there is none).
 *
 * The estimate of a child of a node of bound L, branched upon the candidate
 * b, is that of the best-estimate search of Benichou et al. (1971), i.e., L
 * worsened by
 * \f[
 *  \Psi^d_b \, \delta^d_b + \sum_{v \neq b} \min \{ \Psi^0_v \,
 *  \delta^0_v \, , \, \Psi^1_v \, \delta^1_v \}
 * \f]
 * for the direction d of the child and the distances \f$ \delta^d_v \f$ of
 * the candidates v: the cheapest rounding of every fractional candidate,
 * that of b being forced by the child. The ids of the candidates have to
 * name the same variable at all the nodes. */

class PseudoCosts
{
 public:

 /// the priors of the candidate of the given index, down and up, per unit
 /// of distance (NaN for none)
 using Prior = std::function< std::array< double , 2 >( std::size_t ) >;

 /// record that branching on the candidate id in the direction dir, its
 /// value moving by moved > 0, worsened the bound by worse (taken as 0 if
 /// negative)
 void record( std::size_t id , int dir , double moved , double worse );

 /// the pseudocost of the candidate id in the direction dir, or prior (if
 /// not NaN) if it has no observation in that direction
 [[nodiscard]] double get( std::size_t id , int dir , double prior ) const;

 /// whether the candidate id has an observation in the direction dir
 [[nodiscard]] bool known( std::size_t id , int dir ) const;

 /// the number of the observations recorded so far
 [[nodiscard]] unsigned long observations( void ) const {
  return( f_all[ 0 ].n + f_all[ 1 ].n );
  }

 /// the estimates of the "down" and of the "up" child of the candidate b of
 /// cands, at a node of bound node of a problem to be maximized if maximize;
 /// prior, if given, gives the priors of the candidates by their index
 [[nodiscard]] std::array< double , 2 > estimates(
				const std::vector< BranchCandidate > & cands ,
				std::size_t b , double node , bool maximize ,
				const Prior & prior = nullptr ) const;

 private:

 /// a sum of observations and their number
 struct Sum {
  double s = 0;
  unsigned long n = 0;
  };

 std::unordered_map< std::size_t , std::array< Sum , 2 > > f_sums;  ///< by id
 std::array< Sum , 2 > f_all;                 ///< of all the candidates

 };  // end( class( PseudoCosts ) )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS LearnedPseudoCosts -------------------------*/
/*--------------------------------------------------------------------------*/
/// the priors of the pseudocosts read from a file
/** Two models of LearnedBranching, written one after the other in the same
 * file: the first gives the pseudocost "down" of a candidate out of its
 * features [see BranchingFeatures::features()], the second the "up" one,
 * both over \f$ \max \{ | L | , 1 \} \f$ for the bound L of the node, so
 * that they do not depend on the scale of the data; a negative output counts
 * as 0. They are learned, e.g., from the bounds of the children that strong
 * branching computes [see StrongBranching::write_log()]. */

class LearnedPseudoCosts
{
 public:

 /// read the two models from the given file; throws if they are not there
 void load( const std::string & file );

 /// whether the models have been read
 [[nodiscard]] bool empty( void ) const { return( f_down.empty() ); }

 /// the priors, down and up, of a candidate with the given features at a
 /// node of bound node
 [[nodiscard]] std::array< double , 2 > prior(
		      const BranchingFeatures::Features & x , double node ) const;

 private:

 LearnedBranching f_down;   ///< the model of the pseudocosts "down"
 LearnedBranching f_up;     ///< the model of the pseudocosts "up"

 };  // end( class( LearnedPseudoCosts ) )

/** @} end( group( BranchingRule_CLASSES ) ) */

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* BranchingRule.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File BranchingRule.h --------------------------*/
/*--------------------------------------------------------------------------*/
