include ../../mitool/Makefile.in
include ../../mitool/milib/Makefile.in
include ../../mitool/miblaslib/Makefile.in

INCLUDE = $(MIBLASLIB_INC) $(MILIB_INC)

TARGET := libsrtmathlib.a

OBJ := srtmathlib.o

all: $(TARGET) 

libsrtmathlib.a: ${OBJ}
	$(AR) -r $@ $(OBJ)
	ranlib $@

libsrtmathlib.so: ${OBJ}
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
