/*--------------------------------------------------------------------------*/
/*--------------------------- File GroupChange.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class GroupChange, which implements the
 * Change concept [see Change.h] for a group of Change applied as one: the
 * natural representation of a branching decision that changes more than one
 * thing at a time (say, fixes several variables at once). Not used by
 * BranchAndXSolver yet.
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

#ifndef __GroupChange
 #define __GroupChange
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Change.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS GroupChange -----------------------------*/
/*--------------------------------------------------------------------------*/
/// concrete class for a group of Change applied as one
/** The class GroupChange derives from Change and provides a STL container
 * of pointers to Change, applied in order as a single Change. Since a
 * GroupChange is a Change, nothing prevents some of its sub-Change from
 * being GroupChange themselves, allowing arbitrarily nested Change. The
 * GroupChange owns its sub-Change. */

class GroupChange : public Change {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

 explicit GroupChange() : Change() {}     ///< constructor: does nothing

 /// destructor: deletes the owned sub-Changes

 ~GroupChange() override {
  for( auto chg : v_sub_Changes )
   delete chg;
  }

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

 /// apply all the sub-Changes, in order, to the given Block
 /** Apply all the sub-Changes, in their order, to the given Block. With
  * \p doUndo true, the returned Change is the GroupChange of the undo
  * Changes of the sub-Changes in *reverse* order, so that applying it
  * brings the Block back to the state prior to the call. */

 Change * apply( Block * block , bool doUndo = false ,
                 ModParam issueMod = eNoBlck ,
                 ModParam issueAMod = eNoBlck ) override {
  GroupChange * undo = doUndo ? new GroupChange() : nullptr;
  for( auto chg : v_sub_Changes ) {
   auto u = chg->apply( block , doUndo , issueMod , issueAMod );
   if( doUndo )
    undo->v_sub_Changes.push_front( u );
   }
  return( undo );
  }

/*--------------------------------------------------------------------------*/

 /// add a new Change at the end of the group; ownership is transferred

 void add( Change * chg ) { v_sub_Changes.push_back( chg ); }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

 /// returns a const reference to the list of sub-Changes

 [[nodiscard]] const std::list< Change * > & sub_Changes( void ) const {
  return( v_sub_Changes );
  }

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR LOADING, PRINTING & SAVING THE Change ----------*/
/*--------------------------------------------------------------------------*/

 /// deserialize a GroupChange out of a netCDF::NcGroup: not implemented yet

 void deserialize( const netCDF::NcGroup & group ) override {
  throw( std::logic_error( "GroupChange::deserialize: not implemented yet"
                           ) );
  }

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED METHODS ----------------------------*/
/*--------------------------------------------------------------------------*/

 /// print the GroupChange, i.e., all its sub-Changes

 void print( std::ostream & output ) const override {
  output << "GroupChange:" << std::endl;
  for( const auto chg : v_sub_Changes )
   output << *chg << std::endl;
  }

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS ----------------------------*/
/*--------------------------------------------------------------------------*/

 std::list< Change * > v_sub_Changes;   ///< the std::list of sub-Changes

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 SMSpp_insert_in_factory_h;       // insert GroupChange in the Change factory

/*--------------------------------------------------------------------------*/

 };  // end( class( GroupChange ) )

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* GroupChange.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File GroupChange.h --------------------------*/
/*--------------------------------------------------------------------------*/
