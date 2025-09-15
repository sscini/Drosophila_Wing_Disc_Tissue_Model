-include ../makefile.init

RM := rm -rf

# pull in your per-subdir .cpp / .cu ? .o rules (if you have them)
-include subdir.mk

# All of your .o’s come from subdir.mk + these extras:
OBJS += XMLParser.o \
        gradientRelax.o        # <--- add your new object

# if you have dependency files, keep pulling them in
ifneq ($(MAKECMDGOALS),clean)
  ifneq ($(strip $(C++_DEPS)),)
    -include $(C++_DEPS)
  endif
  ifneq ($(strip $(C_DEPS)),)
    -include $(C_DEPS)
  endif
  ifneq ($(strip $(CC_DEPS)),)
    -include $(CC_DEPS)
  endif
  ifneq ($(strip $(CPP_DEPS)),)
    -include $(CPP_DEPS)
  endif
  ifneq ($(strip $(CXX_DEPS)),)
    -include $(CXX_DEPS)
  endif
  ifneq ($(strip $(C_UPPER_DEPS)),)
    -include $(C_UPPER_DEPS)
  endif
endif

-include ../makefile.defs

#------------------------------------------------------------------------------
# Compiler flags
#------------------------------------------------------------------------------
NVCCFLAGS = -std=c++14

#------------------------------------------------------------------------------
# Implicit rule: compile any .cu ? .o with NVCC (device code)
#------------------------------------------------------------------------------
%.o: %.cu gradientRelax.h
	@echo "Compiling CUDA $< ? $@"
	$(NVCC) $(NVCCFLAGS) -dc -o $@ $<

#------------------------------------------------------------------------------
# All Target
#------------------------------------------------------------------------------
all: tissue-model

#------------------------------------------------------------------------------
# Link
#------------------------------------------------------------------------------
tissue-model: $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: NVCC C++ Linker'
	nvcc $(ILIBS1) $(LIBS) $(OBJS) $(NVCCFLAGS) -o $@
	@echo 'Finished building target: $@'
	@echo

#------------------------------------------------------------------------------
# Clean
#------------------------------------------------------------------------------
clean:
	-$(RM) $(OBJS) tissue-model
	-@echo

#------------------------------------------------------------------------------
# Other targets you had
#------------------------------------------------------------------------------
post-build:
	-@echo

.PHONY: all clean post-build
.SECONDARY: post-build

-include ../makefile.targets



#-include ../makefile.init
#
#RM := rm -rf
#
## All of the sources participating in the build are defined here
#-include subdir.mk
#
## Append the XMLParser object so that XMLParser.cpp gets compiled.
#OBJS += XMLParser.o
#
#ifneq ($(MAKECMDGOALS),clean)
#ifneq ($(strip $(C++_DEPS)),)
#-include $(C++_DEPS)
#endif
#ifneq ($(strip $(C_DEPS)),)
#-include $(C_DEPS)
#endif
#ifneq ($(strip $(CC_DEPS)),)
#-include $(CC_DEPS)
#endif
#ifneq ($(strip $(CPP_DEPS)),)
#-include $(CPP_DEPS)
#endif
#ifneq ($(strip $(CXX_DEPS)),)
#-include $(CXX_DEPS)
#endif
#ifneq ($(strip $(C_UPPER_DEPS)),)
#-include $(C_UPPER_DEPS)
#endif
#endif
#
#-include ../makefile.defs
#
## Add inputs and outputs from these tool invocations to the build variables 
#
## All Target
#all: virus-model
#
##Flags 
##CXXFLAGS=-O3 -std=c++0x -pg -g -c -Wall
#NVCCFLAGS= -std=c++11
##  -g -G -O2
#
## Tool invocations 
#virus-model: $(OBJS)
#	@echo 'Building target: $@'
#	@echo 'Invoking: NVCC C++ Linker'
#	nvcc $(ILIBS1) $(LIBS) $(OBJS) $(NVCCFLAGS) -o "virus-model"  
#	@echo 'Finished building target: $@'
#	@echo ' ' 
#	$(MAKE) --no-print-directory post-build 
#
## Other Targets
#clean:
#	-$(RM) $(OBJS)$(C++_DEPS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS) virus-model
#	-@echo ' '
#
#post-build:
#	#-mkdir --parents ../../dpd-0.0.6/Debug; cd ../../dpd-0.0.6/Debug; cp ../../dpd-model/Debug/dpd-model dpd-model
#	-@echo ' '
#
#.PHONY: all clean dependents
#.SECONDARY: post-build
#
#-include ../makefile.targets
## DO NOT DELETE
