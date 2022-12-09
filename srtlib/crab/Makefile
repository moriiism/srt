include ../../../mitool/Makefile.in
include ../../../mitool/milib/Makefile.in
include ../../../mitool/miblaslib/Makefile.in
include ../../../mitool/mifitslib/Makefile.in
include ../../../mitool/mirootlib/Makefile.in

INCLUDE = $(MIROOTLIB_INC) $(MIFITSLIB_INC) $(MIBLASLIB_INC) \
$(MILIB_INC) $(CFITSIO_INC) $(ROOT_INC)

TARGET := libsrtlib_crab.a

OBJ := rl_crab.o  rl_statval_crab.o


all: $(TARGET) 

libsrtlib_crab.a: ${OBJ}
	$(AR) -r $@ $(OBJ)
	ranlib $@

libsrtlib_crab.so: ${OBJ}
	$(CC) -shared -g -O2 -o $@ ${OBJ}

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
