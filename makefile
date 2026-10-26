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

BAXSLVOBJ = $(BAXSLVSDR)/obj/BranchAndXSolver.o \
	$(BAXSLVSDR)/obj/BranchingRule.o

BAXSLVINC = -I$(BAXSLVSDR)/include

BAXSLVH   = $(BAXSLVSDR)/include/BranchAndXSolver.h \
	$(BAXSLVSDR)/include/BranchingRule.h

# the C++20 module smspp.branchandxsolver - - - - - - - - - - - - - - - - - -
# [see makefile-mod of SMS++]

$(eval $(call SMSppMODrule,BAXSLV,$(BAXSLVSDR),BranchAndXSolver,\
	smspp.branchandxsolver,SMS++))

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(BAXSLVOBJ) $(BAXSLVSDR)/*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(BAXSLVSDR)/obj/BranchAndXSolver.o: \
	$(BAXSLVSDR)/src/BranchAndXSolver.cpp $(BAXSLVH) $(SMS++OBJ)
	$(CC) -c $(BAXSLVSDR)/src/BranchAndXSolver.cpp -o $@ \
	$(BAXSLVINC) $(SMS++INC) $(SW)

$(BAXSLVSDR)/obj/BranchingRule.o: \
	$(BAXSLVSDR)/src/BranchingRule.cpp $(BAXSLVSDR)/include/BranchingRule.h
	$(CC) -c $(BAXSLVSDR)/src/BranchingRule.cpp -o $@ \
	$(BAXSLVINC) $(SW)


############################ End of makefile #################################
