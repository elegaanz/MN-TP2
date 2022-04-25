all: test_gemv test_complexe4 test_complexe3 test_complexe2 test_complexe test_dot test_dot_dyn test_axpy test_asum test_iamin test_iamax test_nrm2 test_gemm 

INC_DIR=../include

LIB_DIR=../lib

LIBST=-lmnblas
LIBDYN=-lmnblasdyn

OPTIONS_COMPIL  =-Wall -O2 -fPIC -I$(INC_DIR) -g

OPTIONS_LINK_STATIC  = -fopenmp -L$(LIB_DIR) $(LIBST) -lm -g
OPTIONS_LINK_DYN     = -fopenmp -L$(LIB_DIR) $(LIBDYN) -lm -g


test_gemv: test_gemv.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_gemv test_gemv.o flop.o $(OPTIONS_LINK_STATIC)

test_gemm: test_gemm.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_gemm test_gemm.o flop.o $(OPTIONS_LINK_STATIC)

test_dot: test_dot.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_dot test_dot.o flop.o $(OPTIONS_LINK_STATIC)

flop.o: flop.c flop.h
	gcc $(OPTIONS_COMPIL) -c flop.c

test_dot_dyn: test_dot.o flop.o
	gcc -o test_dot_dyn flop.o test_dot.o $(OPTIONS_LINK_DYN)

test_asum: test_asum.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_asum test_asum.o flop.o $(OPTIONS_LINK_STATIC)

test_iamin: test_iamin.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_iamin test_iamin.o flop.o $(OPTIONS_LINK_STATIC)

test_iamax: test_iamax.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_iamax test_iamax.o flop.o $(OPTIONS_LINK_STATIC)

test_%.o: test_%.c $(INC_DIR)/mnblas.h
	gcc $(OPTIONS_COMPIL) -c $<

test_axpy: test_axpy.o flop.o $(LIB_DIR)/libmnblas.a
	gcc -o test_axpy test_axpy.o flop.o $(OPTIONS_LINK_STATIC)

test_nrm2: test_nrm2.o flop.o
	gcc -o test_nrm2 test_nrm2.o flop.o $(OPTIONS_LINK_STATIC)

test_complexe.o: test_complexe.c  $(INC_DIR)/complexe.h
	gcc $(OPTIONS_COMPIL) -c test_complexe.c

test_complexe: test_complexe.o flop.o
	gcc -o test_complexe test_complexe.o flop.o $(OPTIONS_LINK_STATIC)

test_complexe2.o: test_complexe2.c  $(INC_DIR)/complexe2.h
	gcc $(OPTIONS_COMPIL) -c test_complexe2.c

test_complexe2: test_complexe2.o flop.o
	gcc -o test_complexe2 test_complexe2.o flop.o $(OPTIONS_LINK_STATIC)


test_complexe3.o: test_complexe3.c  $(INC_DIR)/complexe2.h
	gcc $(OPTIONS_COMPIL) -c test_complexe3.c

test_complexe3: test_complexe3.o flop.o
	gcc -o test_complexe3 test_complexe3.o flop.o $(OPTIONS_LINK_STATIC)

test_complexe4.o: test_complexe4.c  $(INC_DIR)/complexe2.h
	gcc $(OPTIONS_COMPIL) -c test_complexe4.c

test_complexe4: test_complexe4.o flop.o
	gcc -o test_complexe4 test_complexe4.o flop.o $(OPTIONS_LINK_STATIC)


clean:
	rm -f *.o test_dot test_dot_dyn test_complexe test_complexe2 test_complexe3 test_complexe4 test_axpy test_asum test_iamin test_nrm2 test_iamax test_gemm *~
