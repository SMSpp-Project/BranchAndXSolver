##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of BranchAndXSolver                                         #
#                                                                            #
#   Note that $(SMS++INC) is assumed to include any -I directive             #
#   corresponding to external libraries needed by SMS++, at least to the     #
#   extent in which they are needed by the parts of SMS++ used by            #
#   BranchAndXSolver, and that $(BKBkINC) / $(BKBkH) / $(BKBkOBJ) come   #
#   from the makefile of the required BinaryKnapsackBlock module.            #
#                                                                            #
#   Input:  $(CC)       = compiler command                                   #
#           $(SW)       = compiler options                                   #
#           $(SMS++INC) = the -I$( core SMS++ directory )                    #
#           $(SMS++OBJ) = the libSMS++ library itself                        #
#           $(BKBkINC)  = the -I$( BinaryKnapsackBlock directory )           #
#           $(BKBkH)    = the .h files of BinaryKnapsackBlock                #
#           $(BKBkOBJ)  = the BinaryKnapsackBlock object(s)                  #
#           $(BAXSLVSDR)  = the directory where the source is                  #
#                                                                            #
#   Output: $(BAXSLVOBJ)  = the final object(s) / library                      #
#           $(BAXSLVH)    = the .h files to include                            #
#           $(BAXSLVINC)  = the -I$( source directory )                        #
#                                                                            #
#                              Filippo Magi                                  #
#                              Donato Meoli                                  #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

BAXSLVOBJ = $(BAXSLVSDR)/obj/BranchAndXSolver.o \
          $(BAXSLVSDR)/obj/GreedyRelaxationSolver.o

BAXSLVINC = -I$(BAXSLVSDR)/include

BAXSLVH   = $(BAXSLVSDR)/include/BranchAndXSolver.h \
          $(BAXSLVSDR)/include/GreedyRelaxationSolver.h \
          $(BAXSLVSDR)/include/ParallelSolver.h

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(BAXSLVOBJ) $(BAXSLVSDR)/*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(BAXSLVSDR)/obj/BranchAndXSolver.o: \
	$(BAXSLVSDR)/src/BranchAndXSolver.cpp $(BAXSLVH) $(BKBkH) $(SMS++OBJ)
	$(CC) -c $(BAXSLVSDR)/src/BranchAndXSolver.cpp -o $@ \
	$(BAXSLVINC) $(BKBkINC) $(SMS++INC) $(SW)

$(BAXSLVSDR)/obj/GreedyRelaxationSolver.o: \
	$(BAXSLVSDR)/src/GreedyRelaxationSolver.cpp $(BAXSLVH) $(BKBkH) $(SMS++OBJ)
	$(CC) -c $(BAXSLVSDR)/src/GreedyRelaxationSolver.cpp -o $@ \
	$(BAXSLVINC) $(BKBkINC) $(SMS++INC) $(SW)

############################ End of makefile #################################
