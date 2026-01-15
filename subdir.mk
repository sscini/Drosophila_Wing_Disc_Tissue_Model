################################################################################
# Sources / objects list
################################################################################

# CUDA sources (adjust paths if your sources are in ../)
CU_SRCS += \
LinearSprings.cu \
VolumeComp.cu \
VolumeSprings.cu \
NodeAdvance.cu \
System.cu \
Utilities.cu \
StrainTensor.cu \
SystemBuilder.cu \
Storage.cu \
main.cu

OBJS += \
LinearSprings.o \
VolumeComp.o \
VolumeSprings.o \
NodeAdvance.o \
System.o \
Utilities.o \
StrainTensor.o \
SystemBuilder.o \
Storage.o \
main.o
