LIBS = 

BINS = cvmap  dca  diff_logl_map  emgist  emgist_case2  \
	emgist_case2_ls  emgist_case2_lsig  eval  mkresp  \
	pmmn  script  simobs  smrmap

SUBDIR = $(LIBS) $(BINS) script  

ALL = $(LIBS) $(BINS)

all:
	for dir in $(ALL); do \
	(cd $$dir; ${MAKE} all); \
	done
test:
	for dir in $(LIBS); do \
	(cd $$dir; ${MAKE} test); \
	done

clean:
	for dir in $(SUBDIR); do \
	(cd $$dir; ${MAKE} clean); \
	done

cleaner:
	-rm -f *~
	for dir in $(SUBDIR); do \
	(cd $$dir; ${MAKE} cleaner); \
	done
