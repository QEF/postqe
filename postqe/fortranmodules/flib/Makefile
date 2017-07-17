
.SUFFIXES :
.SUFFIXES : .o .c .f .f90



.f90.o:
	gfortran -c -fPIC $<


LIBFILES=\
generate_gg_list.o\
vloc_of_g.o\
erf.o\
simpsn.o\
latgen.o\
errore.o\
struct_fact.o


all: libQE.a
	$(dummy-variable)

libQE.a: $(LIBFILES)
	ar -r $@ $?
clean:
	rm *.o
	touch *.f90

