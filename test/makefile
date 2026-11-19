##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of test (for CSPBlock)                                          #
#                                                                            #
#                                Donato Meoli                                #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# module name, and the name of the battery
NAME = CSPBlock_unit_test
BNAME = CSPBlock_battery_test

# basic directory
DIR = .

# compiler
CC = clang++

# common flags
COMMON_SW = -std=c++20 $(SMSppSWERR) -Wno-enum-compare

# debug switches
SW_DEBUG = -g3 $(SMSppSWDBG) -fno-inline $(COMMON_SW)

# production switches
SW_RELEASE = -O3 $(COMMON_SW) -DNDEBUG

# default target - build in release mode
default: SW = $(SW_RELEASE)
default: $(DIR)/$(NAME) $(DIR)/$(BNAME)

# debug target - build in debug mode
debug: SW = $(SW_DEBUG)
debug: $(DIR)/$(NAME) $(DIR)/$(BNAME)

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -

# SMS++
SMS++SDR = ../../SMS++
include $(SMS++SDR)/lib/makefile-inc

# CSPBlock
CSPBkSDR = ..
include $(CSPBkSDR)/makefile-s

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

# object files
MOBJ =  $(CSPBkOBJ) $(SMS++OBJ)

# libraries
MLIB =  $(CSPBkLIB) $(SMS++LIB)

$(DIR)/$(NAME): $(MOBJ) $(DIR)/test.o
	$(CC) -o $(DIR)/$(NAME) $(DIR)/test.o $(MOBJ) $(MLIB) $(SW)

$(DIR)/$(BNAME): $(MOBJ) $(DIR)/battery.o
	$(CC) -o $(DIR)/$(BNAME) $(DIR)/battery.o $(MOBJ) $(MLIB) $(SW)

# dependencies: every .o from its .cpp + every recursively included .h- - - -

# include directives
MINC =  $(SMS++INC) $(CSPBkINC)

# includes
MH =    $(SMS++H) $(CSPBkH)

# compile command
$(DIR)/test.o: $(DIR)/test.cpp $(MH)
	$(CC) -c $(DIR)/test.cpp -o $@ $(MINC) $(SW)

$(DIR)/battery.o: $(DIR)/battery.cpp $(MH)
	$(CC) -c $(DIR)/battery.cpp -o $@ $(MINC) $(SW)

# clean target
clean::
	rm -f $(DIR)/*.o $(DIR)/*~
	rm -f $(DIR)/$(NAME) $(DIR)/$(BNAME)

# phony targets
.PHONY: default debug clean

############################ End of makefile #################################
