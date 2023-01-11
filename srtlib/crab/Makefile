include ../../../mitool/Makefile.in
include ../../../mitool/milib/Makefile.in
include ../../../mitool/miblaslib/Makefile.in
include ../../../mitool/mifitslib/Makefile.in
include ../../../mitool/mirootlib/Makefile.in

CUDABLAS_INC = -I/usr/local/cuda/include

INCLUDE = $(CUDABLAS_INC) $(MIROOTLIB_INC) $(MIFITSLIB_INC) $(MIBLASLIB_INC) \
	$(MILIB_INC) $(CFITSIO_INC) $(ROOT_INC)

TARGET := libsrtlib_crab.a  libsrtlib_crab_cuda.a

OBJ := rl_crab.o  rl_crab_statval.o \
	rl_crab_smth_pf_zal.o

OBJ_CUDA := rl_crab_cuda.o  rl_crab_statval_cuda.o \
	rl_crab_smth_pf_zal_cuda.o

all: $(TARGET) 

libsrtlib_crab.a: ${OBJ}
	$(AR) -r $@ $(OBJ)
	ranlib $@

libsrtlib_crab.so: ${OBJ}
	$(CC) -shared -g -O2 -o $@ ${OBJ}

libsrtlib_crab_cuda.a: ${OBJ_CUDA}
	$(AR) -r $@ $(OBJ_CUDA)
	ranlib $@

libsrtlib_crab_cuda.so: ${OBJ_CUDA}
	$(CC) -shared -g -O2 -o $@ ${OBJ_CUDA}

.SUFFIXES:
.SUFFIXES: .o .c .f .F .cc .cu

.c.o:
	$(CCN) -c $<

.F.o : 
	$(FCN) -c $<

.f.o:
	$(FCN) -c $<

.cc.o:
	${CXXN} -c $<

.cu.o:
	${NVCCN} -c $<

install: all
	if [ ! -d $(INST_BIN_DIR) ]; then mkdir $(INST_BIN_DIR); fi
	install ${TARGET} ${INST_BIN_DIR}

clean:
	-rm -f *.o

cleaner: clean
	-rm -f ${TARGET} *% *~ core
