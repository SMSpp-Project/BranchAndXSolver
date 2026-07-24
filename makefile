##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of BranchAndXSolver                                             #
#                                                                            #
#   Note that $(SMS++INC) is assumed to include any -I directive             #
#   corresponding to external libraries needed by SMS++, at least to the     #
#   extent in which they are needed by the parts of SMS++ used by            #
#   BranchAndXSolver, i.e., core SMS++.                                      #
#                                                                            #
#   Input:  $(CC)        = compiler command                                  #
#           $(SW)        = compiler options                                  #
#           $(SMS++INC)  = the -I$( core SMS++ directory )                   #
#           $(SMS++OBJ)  = the libSMS++ library itself                       #
#           $(BAXSLVSDR) = the directory where the source is                 #
#                                                                            #
#   Output: $(BAXSLVOBJ)  = the final object(s) / library                    #
#           $(BAXSLVH)    = the .h files to include                          #
#           $(BAXSLVINC)  = the -I$( source directory )                      #
#                                                                            #
#                              Donato Meoli                                  #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

BAXSLVOBJ = $(BAXSLVSDR)/obj/BranchAndXSolver.o

BAXSLVINC = -I$(BAXSLVSDR)/include

BAXSLVH   = $(BAXSLVSDR)/include/BranchAndXSolver.h

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(BAXSLVOBJ) $(BAXSLVSDR)/*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(BAXSLVSDR)/obj/BranchAndXSolver.o: \
	$(BAXSLVSDR)/src/BranchAndXSolver.cpp $(BAXSLVH) $(SMS++OBJ)
	$(CC) -c $(BAXSLVSDR)/src/BranchAndXSolver.cpp -o $@ \
	$(BAXSLVINC) $(FF_INC) $(SMS++INC) $(SW)


############################ End of makefile #################################
