
QE_DIR = `pwd`/build/q-e

.SUFFIXES :
.SUFFIXES : .o .c .f .f90

.f90.o:
	gfortran -c -fPIC $<

LIBFILES=\
error_handler.o\
generate_gg_list.o\
struct_fact.o

libpostqe: $(LIBFILES)
	cp -p $(QE_DIR)/UtilXlib/libutil.a postqe2.a
	ar -r $@.a $?

pyqe: wrappers_utils.f90 wrappers_funct.f90 wrappers_vloc.f90
	f2py3 -c wrappers_utils.f90 wrappers_funct.f90 wrappers_vloc.f90 -m pyqe \
	-I$(QE_DIR)/Modules/	$(QE_DIR)/clib/c_mkdir.o \
	$(QE_DIR)/clib/copy.o \
	$(QE_DIR)/clib/cptimer.o \
	$(QE_DIR)/clib/customize_signals.o \
	$(QE_DIR)/clib/eval_infix.o \
	$(QE_DIR)/clib/md5_from_file.o \
	$(QE_DIR)/clib/md5.o \
	$(QE_DIR)/clib/memstat.o \
	$(QE_DIR)/clib/ptrace.o \
	$(QE_DIR)/clib/qmmm_aux.o \
	$(QE_DIR)/clib/qsort.o \
	$(QE_DIR)/clib/sockets.o \
	$(QE_DIR)/clib/stack.o \
	-L$(QE_DIR)/Modules/ -lqemod \
	-L$(QE_DIR)/FFTXlib/ -lqefft \
	$(QE_DIR)/PW/src/libpw.a \
	-L$(QE_DIR)/UtilXlib/ -lutil \
        -llapack -lblas \
	libpostqe.a

all: libpostqe pyqe move

move :
	mv *.so ../

clean:
	rm *.o
	touch *.f90
	rm *.so

touch-dummy:
	$(dummy-variable)
