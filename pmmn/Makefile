include ../../mitool/Makefile.in
include ../../mitool/milib/Makefile.in
include ../../mitool/miblaslib/Makefile.in
include ../../mitool/mifitslib/Makefile.in
include ../../mitool/mirootlib/Makefile.in

INCLUDE = $(MIROOTLIB_INC) $(MIFITSLIB_INC) $(MIBLASLIB_INC) $(MILIB_INC) $(CFITSIO_INC) $(ROOT_INC) $(MT_INC)
LDLIBS  = $(MIROOTLIB_LIB) $(MIFITSLIB_LIB) $(MIBLASLIB_LIB) $(MILIB_LIB) $(CFITSIO_LIB) \
	$(ROOT_ALLLIBS) $(LAPACK_LIB) $(BLAS_LIB) $(MT_LIB) $(LIB) -lgfortran

INCLUDE_OPENBLAS = $(MIROOTLIB_INC) $(MIFITSLIB_INC) $(MIBLASLIB_INC) $(MILIB_INC) $(CFITSIO_INC) $(ROOT_INC) $(MT_INC)
LDLIBS_OPENBLAS  = $(MIROOTLIB_LIB) $(MIFITSLIB_LIB) $(MIBLASLIB_LIB) $(MILIB_LIB) $(CFITSIO_LIB) \
	$(ROOT_ALLLIBS) $(LAPACK_LIB) $(OPENBLAS_LIB) $(MT_LIB) $(LIB) -lgfortran


TARGET = pmmn_case2_openblas  pmmn_case2

all: $(TARGET) 

pmmn_case2: pmmn_case2.o  arg_pmmn_case2.o  sub_pmmn_case2.o
	$(CXXN) $@.o arg_$@.o sub_$@.o ${LDLIBS} -o $@

pmmn_case2_openblas: pmmn_case2.o  arg_pmmn_case2.o  sub_pmmn_case2.o
	$(CXXN_OPENBLAS) pmmn_case2.o arg_pmmn_case2.o sub_pmmn_case2.o ${LDLIBS_OPENBLAS} -o $@

pmmn_case2_test_openblas: pmmn_case2.o  arg_pmmn_case2.o  sub_pmmn_case2.o
	$(CXXN_OPENBLAS) pmmn_case2.o arg_pmmn_case2.o sub_pmmn_case2.o ${LDLIBS_OPENBLAS} -o $@


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
