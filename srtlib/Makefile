include ../../mitool/Makefile.in
include ../../mitool/milib/Makefile.in
include ../../mitool/miblaslib/Makefile.in
include ../../mitool/mifitslib/Makefile.in
# include ../../mitool/mirootlib/Makefile.in

CUDABLAS_INC = -I/usr/local/cuda/include

INCLUDE = $(CUBLAS_INC) $(MIFITSLIB_INC) $(MIBLASLIB_INC) \
	$(MILIB_INC) $(CFITSIO_INC)

# $(MIROOTLIB_INC) 
# $(ROOT_INC)

TARGET := libsrtlib.a  libsrtlib_cuda.a

OBJ := smth_zal.o  rl.o  rl_statval.o  rl_det2.o

OBJ_CUDA := smth_zal.o  smth_zal_cuda.o  rl.o  rl_statval.o  rl_det2.o

all: $(TARGET) 

libsrtlib.a: ${OBJ}
	$(AR) -r $@ $(OBJ)
	ranlib $@

libsrtlib.so: ${OBJ}
	$(CC) -shared -g -O2 -o $@ ${OBJ}

libsrtlib_cuda.a: ${OBJ_CUDA}
	$(AR) -r $@ $(OBJ_CUDA)
	ranlib $@

libsrtlib_cuda.so: ${OBJ_CUDA}
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
