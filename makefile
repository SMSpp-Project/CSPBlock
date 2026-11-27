##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of CSPBlock                                                     #
#                                                                            #
#   Note that $(SMS++INC) is assumed to include any -I directive             #
#   corresponding to external libraries needed by SMS++, at least to the     #
#   extent in which they are needed by the parts of SMS++ used by            #
#   CSPBlock.                                                                #
#                                                                            #
#   Input:  $(CC)        = compiler command                                  #
#           $(SW)        = compiler options                                  #
#           $(SMS++INC)  = the -I$( core SMS++ directory )                   #
#           $(SMS++OBJ)  = the libSMS++ library itself                       #
#           $(CSPBkSDR)  = the directory where the source is                 #
#                                                                            #
#   Output: $(CSPBkOBJ)  = the final object(s) / library                     #
#           $(CSPBkH)    = the .h files to include                           #
#           $(CSPBkINC)  = the -I$( source directory )                       #
#                                                                            #
#                             Antonio Frangioni                              #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

CSPBkOBJ = $(CSPBkSDR)/obj/CSPBlock.o \
           $(CSPBkSDR)/obj/FlatZincCSPSolver.o

CSPBkINC = -I$(CSPBkSDR)/include

CSPBkH   = $(CSPBkSDR)/include/CSPBlock.h \
           $(CSPBkSDR)/include/CSPGlobalConstraints.h \
           $(CSPBkSDR)/include/FlatZincCSPSolver.h

# the C++20 module [see makefile-mod of SMS++] - - - - - - - - - - - - - - -
# the prefixes of the modules of the dependencies (e.g., UCBck) go in place
# of SMS++

$(eval $(call SMSppMODrule,CSPBk,$(CSPBkSDR),CSPBlock,\
	smspp.cspblock,SMS++))

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(CSPBkOBJ) $(CSPBkSDR)/*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(CSPBkSDR)/obj/CSPBlock.o: $(CSPBkSDR)/src/CSPBlock.cpp \
	$(CSPBkSDR)/include/CSPBlock.h \
	$(CSPBkSDR)/include/CSPGlobalConstraints.h $(SMS++H) $(SMS++OBJ)
	$(CC) -c $(CSPBkSDR)/src/CSPBlock.cpp -o $@ \
	$(CSPBkINC) $(SMS++INC) $(SW)

$(CSPBkSDR)/obj/FlatZincCSPSolver.o: \
	$(CSPBkSDR)/src/FlatZincCSPSolver.cpp $(CSPBkH) $(SMS++OBJ)
	$(CC) -c $(CSPBkSDR)/src/FlatZincCSPSolver.cpp -o $@ \
	$(CSPBkINC) $(SMS++INC) $(SW)

########################## End of makefile ###################################
