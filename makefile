-include ../makefile.init

RM := rm -rf

# Pull in subdir.mk which defines OBJS and rules (we will fix subdir.mk too)
-include subdir.mk

# Add extra .o not listed in subdir.mk (if these are real sources in this dir)
OBJS += XMLParser.o gradientRelax.o

# ----------------------------
# Toolchain
# ----------------------------
CXX  ?= g++
NVCC ?= nvcc

# ----------------------------
# Includes / libs
# ----------------------------
current_dir := $(shell pwd)

CUDA_INC := -I/opt/linux/rocky/8.x/x86_64/pkgs/cuda/12.1/include/
PUGI_LIB := -L$(current_dir)/pugixml/lib64 -lpugixml

# If you have other libs, append here:
LIBS  += $(PUGI_LIB)
ILIBS += $(CUDA_INC)

# ----------------------------
# Flags
# ----------------------------
CXXFLAGS  += -O3 -std=c++14 -Wall -Wextra -pthread
NVCCFLAGS += -O3 -std=c++14

# IMPORTANT:
# enable double atomicAdd by compiling for sm_60+ (Pascal+).
# Use a fatbinary so it runs on multiple nodes (safe on clusters).
NVCCFLAGS += \
  -gencode arch=compute_60,code=sm_60 \
  -gencode arch=compute_70,code=sm_70 \
  -gencode arch=compute_80,code=sm_80 \
  -gencode arch=compute_86,code=sm_86 \
  -gencode arch=compute_86,code=compute_86

# (Optional but recommended) show warnings in device code
NVCCFLAGS += -Xcompiler="-Wall -Wextra"

# ----------------------------
# Default target
# ----------------------------
all: tissue-model

# ----------------------------
# C++ compile rule
# ----------------------------
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(ILIBS) -c $< -o $@

# ----------------------------
# CUDA compile rule
# ----------------------------
%.o: %.cu
	$(NVCC) $(NVCCFLAGS) $(ILIBS) -dc $< -o $@

# ----------------------------
# Link
# ----------------------------
tissue-model: $(OBJS)
	@echo 'Building target: $@'
	$(NVCC) $(NVCCFLAGS) $(OBJS) $(LIBS) -o $@
	@echo 'Finished building target: $@'

# ----------------------------
# Clean
# ----------------------------
clean:
	-$(RM) $(OBJS) tissue-model *.d
	-@echo

.PHONY: all clean
-include ../makefile.targets
