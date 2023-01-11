include ../../../mitool/Makefile.in
include ../../../mitool/milib/Makefile.in
include ../../../mitool/miblaslib/Makefile.in
include ../../../mitool/mifitslib/Makefile.in
include ../../../mitool/mirootlib/Makefile.in
include ../../Makefile.in

CUDABLAS_INC = -I/usr/local/cuda/include
CUDABLAS_LIB = -L/usr/local/cuda/lib64 -lcublas -lcudart

INCLUDE = ${SRTLIB_CRAB_INC} ${SRTLIB_INC} \
	$(MIROOTLIB_INC) $(MIFITSLIB_INC) $(MIBLASLIB_INC) \
	$(MILIB_INC) $(CFITSIO_INC) $(ROOT_INC) $(MT_INC)
LDLIBS  = ${SRTLIB_CRAB_LIB} ${SRTLIB_LIB} \
	$(MIROOTLIB_LIB) $(MIFITSLIB_LIB) $(MIBLASLIB_LIB) \
	$(MILIB_LIB) $(CFITSIO_LIB) $(ROOT_ALLLIBS) $(LAPACK_LIB) $(BLAS_LIB) \
	$(MT_LIB) $(LIB) -lgfortran

INCLUDE_OPENBLAS = ${SRTLIB_CRAB_INC} ${SRTLIB_INC} \
	$(MIROOTLIB_INC) $(MIFITSLIB_INC) \
	$(MIBLASLIB_INC) $(MILIB_INC) $(CFITSIO_INC) $(ROOT_INC) $(MT_INC)
LDLIBS_OPENBLAS  = ${SRTLIB_CRAB_LIB} ${SRTLIB_LIB} \
	$(MIROOTLIB_LIB) $(MIFITSLIB_LIB) \
	$(MIBLASLIB_LIB) $(MILIB_LIB) $(CFITSIO_LIB) $(ROOT_ALLLIBS) \
	$(LAPACK_LIB) $(OPENBLAS_LIB) $(MT_LIB) $(LIB) -lgfortran

INCLUDE_CUDA = ${SRTLIB_CRAB_INC} ${SRTLIB_INC} \
	$(CUDABLAS_INC) \
	$(MIROOTLIB_INC) $(MIFITSLIB_INC) $(MIBLASLIB_INC) \
	$(MILIB_INC) $(CFITSIO_INC) $(ROOT_INC) $(MT_INC)
LDLIBS_CUDA  = ${SRTLIB_CRAB_LIB} ${SRTLIB_LIB} \
	$(CUDABLAS_LIB) \
	$(MIROOTLIB_LIB) $(MIFITSLIB_LIB) $(MIBLASLIB_LIB) \
	$(MILIB_LIB) $(CFITSIO_LIB) $(ROOT_ALLLIBS) $(LAPACK_LIB) $(BLAS_LIB) \
	$(MT_LIB) $(LIB) -lgfortran



TARGET = richlucy_smth_pf_zal

all: $(TARGET) 

richlucy_smth_pf_zal: richlucy_smth_pf_zal.o  arg_richlucy_smth_pf_zal.o
	$(CXXN) $@.o arg_$@.o ${LDLIBS_OPENBLAS} -o $@

.SUFFIXES:
.SUFFIXES: .o .c .f .F .cc

.c.o:
	$(CCN) -c $<

.F.o : 
	$(FCN) -c $<

.f.o:
	$(FCN) -c $<

.cc.o:
	${CXXN} -c $<

install: all
	if [ ! -d $(INST_BIN_DIR) ]; then mkdir $(INST_BIN_DIR); fi
	install ${TARGET} ${INST_BIN_DIR}

clean:
	-rm -f *.o

cleaner: clean
	-rm -f ${TARGET} *% *~ core
