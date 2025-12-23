################################################################################
# Automatically-generated file. Do not edit!
################################################################################
CXX := gcc
NVCC := nvcc
CFLAGS := -pthread -std=c++11 -Wall -Wextra
MFLAGS := -f
NVCCFLAGS := -std=c++11

current_dir := $(shell pwd)
LIBS:=  -lpugixml -L/$(current_dir)/pugixml/lib64
#-lgsl -lgslcblas

ILIBS_cuda8 = -I/opt/linux/rocky/8.x/x86_64/pkgs/cuda/12.1/include/
ILIBS_cuda9 := -I/opt/linux/centos/7.x/x86_64/pkgs/cuda/9.1/include/
ILIBS_cuda11_2 := -I/opt/linux/rocky/8.x/x86_64/pkgs/cuda/12.1/include/

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../LinearSprings.cu \
../VolumeComp.cu \
../VolumeSprings.cu \
../NodeAdvance.cu \
../System.cu \
../Utilities.cu \
../StrainTensor.cu \
../SystemBuilder.cu \
../Storage.cu \
../main.cu


# this is a variable
OBJS += \
./LinearSprings.o \
./VolumeComp.o \
./VolumeSprings.o \
./NodeAdvance.o \
./System.o \
./Utilities.o \
./StrainTensor.o \
./SystemBuilder.o \
./Storage.o \
./main.o


CPP_DEPS += \
./LinearSprings.d \
./VolumeComp.d \
./VolumeSprings.d \
./NodeAdvance.d \
./System.d \
./Utilities.d \
./StrainTensor.d \
./SystemBuilder.d \
./Storage.d \
./main.d

#need o have ILIBS2
#cpp files
%.o : ./%.cpp 
	 $(CXX) $(CFLAGS) $(ILIBS_cuda8) $(LIBS) -o $@ -c $^ 

	
#cuda files
%.o : ./%.cu 
	$(NVCC) $(NVCCFLAGS) $(ILIBS_cuda11_2) -dc -o $@ $^
