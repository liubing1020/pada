MODEL ?= egfngf
METHOD ?= centralized
ARCH ?= sm_20
REAL ?= double
RKmethod ?= RK4
TOP ?= binningOnGPU

CUDA_PATH=/usr/local/cuda

all: bin/${MODEL}_${METHOD}_${ARCH}_${REAL}_${RKmethod}/${MODEL}

include ${METHOD}/Makefile

bin/${MODEL}_${METHOD}_${ARCH}_${REAL}_${RKmethod}/tmp:
	mkdir -p bin/${MODEL}_${METHOD}_${ARCH}_${REAL}_${RKmethod}/tmp

bin/${MODEL}_${METHOD}_${ARCH}_${REAL}_${RKmethod}/${MODEL}: bin/${MODEL}_${METHOD}_${ARCH}_${REAL}_${RKmethod}/tmp rk/${RKmethod}/rk.cu ${TOP}/top.cu bin/${MODEL}_${METHOD}_${ARCH}_${REAL}_${RKmethod}/simulation.cu models/${MODEL}.model config/${MODEL}_${ARCH}_${REAL}.${METHOD} ${TOP}/WriteFiles.c ${TOP}/incBins.cu ${TOP}/binning_structure.h ${TOP}/config_multigpu.h rng/GenRandom.cu
	cp -f  models/${MODEL}.model bin/${MODEL}_${METHOD}_${ARCH}_${REAL}_${RKmethod}/model.cu
	cp -f config/${MODEL}_${ARCH}_${REAL}.${METHOD} bin/${MODEL}_${METHOD}_${ARCH}_${REAL}_${RKmethod}/config.cu
	nvcc -I . -Dreal=${REAL} -I ${TOP} -I bin/${MODEL}_${METHOD}_${ARCH}_${REAL}_${RKmethod} -I ${CUDA_PATH}/include --ptxas-options=-v ${TOP}/top.cu -I rk/${RKmethod} -o  bin/${MODEL}_${METHOD}_${ARCH}_${REAL}_${RKmethod}/${MODEL} -arch=${ARCH}

.PHONY: clean

clean:
	rm -rf bin

