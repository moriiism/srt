include ../../mitool/Makefile.in
include ../../mitool/milib/Makefile.in
include ../../mitool/miblaslib/Makefile.in
include ../../mitool/mifitslib/Makefile.in
include ../../mitool/mirootlib/Makefile.in

INCLUDE = $(MIROOTLIB_INC) $(MIFITSLIB_INC) $(MIBLASLIB_INC) \
$(MILIB_INC) $(CFITSIO_INC) $(ROOT_INC)

TARGET := libsrtlib_ps.a

OBJ := rl.o  rl_statval.o  sim.o  rl_bg.o  smth.o \
	rl_bg2_smth_em.o  rl_bg2_smth_pm.o  rl_bg2_smth_newton.o

#    
#	

#	fpsrc_smth_bg_mm_em.o  fpsrc_smth_bg_mm_pm.o \
#	fpsrc_smth_bg_mm_newton.o

#	fpsrc_smth_bg_newton.o  \
#	fpsrc_smth_bg_statval.o \
#	fpsrc_smth_bg_em.o  fpsrc_smth_bg_dc.o  fpsrc_smth_bg_pm.o  \
#	fpsrc_smth_bg_pmproj_em.o  fpsrc_smth_bg_pmproj_dc.o \
#	fpsrc_smth_bg_pmproj_pm.o \
#	

#sub_smooth.o  sim.o \
#	fpsrc_smth_bg_newton.o  \
#	fpsrc_smth_bg_statval.o \
#	fpsrc_smth_bg_em.o  fpsrc_smth_bg_dc.o  fpsrc_smth_bg_pm.o  \
#	fpsrc_smth_bg_newton_phi0.o 

# sub.o  sub_pm.o  sub_newton.o  
#	fpsrc_smth_bg_pmproj_em.o  fpsrc_smth_bg_pmproj_dc.o \
#	fpsrc_smth_bg_pmproj_pm.o


all: $(TARGET) 

libsrtlib_ps.a: ${OBJ}
	$(AR) -r $@ $(OBJ)
	ranlib $@

libsrtlib_ps.so: ${OBJ}
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
