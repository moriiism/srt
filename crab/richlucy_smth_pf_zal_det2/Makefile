include ../../../mitool/Makefile.in
include ../../../mitool/milib/Makefile.in
include ../../../mitool/miblaslib/Makefile.in
include ../../../mitool/mifitslib/Makefile.in
include ../../../mitool/mirootlib/Makefile.in
include ../../Makefile.in

CUDABLAS_INC = -I/usr/local/cuda/include
CUDABLAS_LIB = -L/usr/local/cuda/lib64 -lcublas -lcudart

INCLUDE = ${SRTLIB_CRAB_INC} ${SRTLIB_INC} ${SRT_MATH_LIB_INC} \
	$(MIFITSLIB_INC) $(MIBLASLIB_INC) \
	$(MILIB_INC) $(CFITSIO_INC)
LDLIBS  = ${SRTLIB_CRAB_LIB} ${SRTLIB_LIB} ${SRT_MATH_LIB_LIB} \
	$(MIFITSLIB_LIB) $(MIBLASLIB_LIB) \
	$(MILIB_LIB) $(CFITSIO_LIB) $(LAPACK_LIB) $(BLAS_LIB) \
	$(LIB) -lgfortran

INCLUDE_OPENBLAS = ${SRTLIB_CRAB_INC} ${SRTLIB_INC} ${SRT_MATH_LIB_INC} \
	$(MIFITSLIB_INC) \
	$(MIBLASLIB_INC) $(MILIB_INC) $(CFITSIO_INC)
LDLIBS_OPENBLAS  = ${SRTLIB_CRAB_LIB} ${SRTLIB_LIB} ${SRT_MATH_LIB_LIB} \
	$(MIFITSLIB_LIB) \
	$(MIBLASLIB_LIB) $(MILIB_LIB) $(CFITSIO_LIB) \
	$(LAPACK_LIB) $(OPENBLAS_LIB) $(LIB) -lgfortran

INCLUDE_CUDA = ${SRTLIB_CRAB_INC} ${SRTLIB_INC} ${SRT_MATH_LIB_INC} \
	$(CUDABLAS_INC) \
	$(MIFITSLIB_INC) $(MIBLASLIB_INC) \
	$(MILIB_INC) $(CFITSIO_INC)
LDLIBS_CUDA  = ${SRTLIB_CRAB_CUDA_LIB} ${SRTLIB_CRAB_LIB} ${SRTLIB_CUDA_LIB} \
	${SRT_MATH_LIB_LIB} $(CUDABLAS_LIB) \
	$(MIFITSLIB_LIB) $(MIBLASLIB_LIB) \
	$(MILIB_LIB) $(CFITSIO_LIB) $(LAPACK_LIB) $(BLAS_LIB) \
	$(LIB) -lgfortran


TARGET = richlucy_smth_pf_zal_det2
#richlucy_smth_pf_zal_cuda  calc_flux0

all: $(TARGET) 

richlucy_smth_pf_zal_det2: richlucy_smth_pf_zal_det2.o  arg_richlucy_smth_pf_zal_det2.o
	$(CXXN) $@.o arg_$@.o ${LDLIBS_OPENBLAS} -o $@

richlucy_smth_pf_zal_cuda: richlucy_smth_pf_zal_cuda.o  arg_richlucy_smth_pf_zal.o
#	$(NVCCN) richlucy_smth_pf_zal_cuda.o  arg_richlucy_smth_pf_zal.o ${LDLIBS_CUDA} -o $@
	$(CXXN) richlucy_smth_pf_zal_cuda.o  arg_richlucy_smth_pf_zal.o ${LDLIBS_CUDA} -o $@

eval_val: eval_val.o  arg_eval_val.o
	$(CXXN) $@.o arg_$@.o ${LDLIBS_OPENBLAS} -o $@

calc_flux0: calc_flux0.o  arg_calc_flux0.o
	$(CXXN) $@.o arg_$@.o ${LDLIBS_OPENBLAS} -o $@


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
