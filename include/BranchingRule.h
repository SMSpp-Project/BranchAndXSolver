/*--------------------------------------------------------------------------*/
/*------------------------- File BranchingRule.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the components that a RelaxationSolver [see
 * ChangeSolver.h] can use to choose the variable to branch upon: the
 * description of a candidate and its features (BranchingFeatures), strong
 * branching with the product score (StrongBranching), and a scoring function
 * of the candidates read from a file (LearnedBranching). They live beside the
 * BranchAndXSolver but know nothing of it: each RelaxationSolver describes
 * its own candidates, decides which component to use (say, by a parameter of
 * its own), and gives StrongBranching the way to evaluate the children of a
 * candidate, so that the same rule serves any relaxation without the
 * RelaxationSolver exposing anything of its representation.
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

#include <ostream>

#include <string>

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

/** @} end( group( BranchingRule_CLASSES ) ) */

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* BranchingRule.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File BranchingRule.h --------------------------*/
/*--------------------------------------------------------------------------*/
