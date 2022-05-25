##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of BranchAndXSolver                                             #
#                                                                            #
#   The makefile takes in input the -I directives for all the external       #
#   libraries needed by BranchAndXSolver, i.e., core SMS++.                  #
#                                                                            #
#   Note that, conversely, $(SMS++INC) is also assumed to include any        #
#   -I directive corresponding to external libraries needed by SMS++, at     #
#   least to the extent in which they are needed by the parts of SMS++       #
#   used by BranchAndXSolver.                                                #
#                                                                            #
#   Input:  $(CC)          = compiler command                                #
#           $(SW)          = compiler options                                #
#           $(SMS++INC)    = the -I$( core SMS++ include directory )         #
#           $(SMS++OBJ)    = the core SMS++ library                          #
#           $(BAXSLVSDR)   = the directory where the source is               #
#                                                                            #
#   Output: $(BAXSLVOBJ)   = the final object(s) / library                   #
#           $(BAXSLVH)     = the .h files to include                         #
#           $(BAXSLVINC)   = the -I$( source directory )                     #
#                                                                            #
#                                VERSION 0.0                                 #
#                               25 - 05 - 2022                               #
#                                                                            #
#                              Antonio Frangioni                             #
#                            Federica Di Pasquale                            #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

BAXSLVOBJ = $(BAXSLVSDR)BranchAndXSolver.o 

BAXSLVINC = -I$(BAXSLVSDR)

BAXSLVH   = $(BAXSLVSDR)BranchAndXSolver.h

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(BAXSLVOBJ) $(BAXSLVSDR)*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(BAXSLVSDR)BranchAndXSolver.o: $(BAXSLVSDR)BranchAndXSolver.cpp \
	$(BAXSLVSDR)BranchAndXSolver.h $(SMS++OBJ) 
	$(CC) -c $*.cpp -o $@ $(BAXSLVINC) $(SMS++INC) $(SW)

########################## End of makefile ###################################
