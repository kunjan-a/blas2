ARCH=64
blas.out : SSYMV.o SGEMV.o STRMV.o SGBMV.o INIT_VARS.o blas2.o
	g++ -m$(ARCH) SSYMV.o SGEMV.o STRMV.o SGBMV.o INIT_VARS.o blas2.o -o blas.out -lrt -ltbb

INIT_VARS.o : INIT_VARS.cpp blas2.h
	g++ -c -g -m$(ARCH) INIT_VARS.cpp -o INIT_VARS.o

SSYMV.o : SSYMV.cpp blas2.h
	g++ -c -g -m$(ARCH) SSYMV.cpp -o SSYMV.o

SGEMV.o : SGEMV.cpp blas2.h
	g++ -c -g -m$(ARCH) SGEMV.cpp -o SGEMV.o

STRMV.o : STRMV.cpp blas2.h
	g++ -c -g -m$(ARCH) STRMV.cpp -o STRMV.o

SGBMV.o : SGBMV.cpp blas2.h
	g++ -c -g -m$(ARCH) SGBMV.cpp -o SGBMV.o

blas2.o : blas2.cpp blas2.h
	g++ -c -g -m$(ARCH) blas2.cpp -o blas2.o

TBB : SGER_tbb.o SSYMV_tbb.o SGEMV_tbb.o STRMV_tbb.o SGBMV_tbb.o INIT_VARS.o blas2.o
	g++ -m$(ARCH) SGER.o SSYMV.o SGEMV.o STRMV.o SGBMV.o INIT_VARS.o blas2.o -o blas.out -lrt -ltbb

SGER_tbb.o : SGER_tbb.cpp blas2.h
	g++ -c -g -m$(ARCH) SGER_tbb.cpp -o SGER.o

SSYMV_tbb.o : SSYMV_tbb.cpp blas2.h
	g++ -c -g -m$(ARCH) SSYMV_tbb.cpp -o SSYMV.o

SGEMV_tbb.o : SGEMV_tbb.cpp blas2.h
	g++ -c -g -m$(ARCH) SGEMV_tbb.cpp -o SGEMV.o

STRMV_tbb.o : STRMV.cpp blas2.h
	g++ -c -g -m$(ARCH) STRMV_tbb.cpp -o STRMV.o

SGBMV_tbb.o : SGBMV.cpp blas2.h
	g++ -c -g -m$(ARCH) SGBMV.cpp -o SGBMV.o


MKL : SSYMV_tbb.o SGEMV_tbb.o STRMV_mkl.o SGBMV_mkl.o INIT_VARS.o blas2.o
	g++ -m$(ARCH) SSYMV.o SGEMV.o STRMV.o SGBMV.o INIT_VARS.o blas2.o -o blas.out -lrt -ltbb -L$MKLROOT/lib -Wl,--start-group -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -Wl,--end-group -fopenmp -lpthread

SSYMV_mkl.o : SSYMV.cpp blas2.h
	g++ -c -g -m$(ARCH) SSYMV.cpp -o SSYMV.o

SGEMV_mkl.o : SGEMV_mkl.cpp blas2.h
	g++ -c -g -m$(ARCH) SGEMV_mkl.cpp -o SGEMV.o

STRMV_mkl.o : STRMV.cpp blas2.h
	g++ -c -g -m$(ARCH) STRMV.cpp -o STRMV.o

SGBMV_mkl.o : SGBMV.cpp blas2.h
	g++ -c -g -m$(ARCH) SGBMV.cpp -o SGBMV.o
clean :
	$(RM) *.o *.out

