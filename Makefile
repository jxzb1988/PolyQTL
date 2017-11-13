CC=g++
DIC=$(PWD)
CFLAGS=-c -Wall -g -O3  -I $(DIC)  
LDFLAGS= -DMKL_ILP64 -m64  -fopenmp  -DARMA_DONT_USE_WRAPPER -llapack -lblas -lgslcblas  -lgsl    -DWITH_LAPACK  -std=c++11  -I  $(DIC)/armadillo-8.100.1/include  -I/usr/local/include  -lz 

#LDFLAGS= -DMKL_ILP64 -m64  -fopenmp  -DARMA_DONT_USE_WRAPPER -llapack -lblas -lgslcblas  -lgsl    -DWITH_LAPACK  -std=c++0x  -I  $(DIC)/armadillo-8.100.1/include  -I/usr/local/include  -lz

 
SOURCES1=conditional_function.cpp MeQTLPolyG.cpp PostCal.cpp Util.cpp TopKSNP.cpp   gemma.cpp  gemma_gzstream.cpp gemma_io.cpp gemma_lapack.cpp gemma_lmm.cpp gemma_mathfunc.cpp gemma_param.cpp 

EXECUTABLE1=PolyQTL
	
$(EXECUTABLE1): $(SOURCES1) 
	$(CC) $(SOURCES1)   $(LDFLAGS) -o $@
clean:
	rm PolyQTL
