include ../../mitool/Makefile.in
include ../../mitool/milib/Makefile.in
include ../../mitool/miblaslib/Makefile.in
include ../../mitool/mifitslib/Makefile.in
include ../../mitool/mirootlib/Makefile.in
include ../Makefile.in

INCLUDE = ${SRTLIB_INC} ${SRT_MATH_LIB_INC} \
	$(MIFITSLIB_INC) \
	$(MIBLASLIB_INC) $(MILIB_INC) $(CFITSIO_INC)
LDLIBS = ${SRTLIB_LIB} ${SRT_MATH_LIB_LIB} \
	$(MIFITSLIB_LIB) \
	$(MIBLASLIB_LIB) $(MILIB_LIB) $(CFITSIO_LIB) \
	$(LAPACK_LIB) $(BLAS_LIB) \
	$(LIB) -lgfortran

INCLUDE_OPENBLAS = ${SRTLIB_INC} ${SRT_MATH_LIB_INC} \
	$(MIFITSLIB_INC) \
	$(MIBLASLIB_INC) $(MILIB_INC) $(CFITSIO_INC)
LDLIBS_OPENBLAS = ${SRTLIB_LIB} ${SRT_MATH_LIB_LIB} \
	$(MIFITSLIB_LIB) \
	$(MIBLASLIB_LIB) $(MILIB_LIB) $(CFITSIO_LIB) \
	$(LAPACK_LIB) $(OPENBLAS_LIB) \
	$(LIB) -lgfortran

TARGET = richlucy

all: $(TARGET) 

richlucy: richlucy.o  arg_richlucy.o
	$(CXXN) $@.o arg_$@.o ${LDLIBS} -o $@

richlucy_openblas: richlucy.o  arg_richlucy.o
	$(CXXN) richlucy.o arg_richlucy.o ${LDLIBS_OPENBLAS} -o $@


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
