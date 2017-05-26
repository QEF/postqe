
.SUFFIXES :
.SUFFIXES : .o .c .f .f90



.f90.o:
	gfortran -c -fPIC $<


libQE.a: touch-dummy
	(cd flib ; make all)

pyQ:  wrappers_and_drivers.f90 libQE.a
	f2py -c wrappers_and_drivers.f90 -m pyQ flib/libQE.a  
	

funct: funct.f90
	f2py -c -m funct funct.f90 functionals.f90 lsda_functionals.f90 expint.f90 erf.f90
all: pyQ funct

move :
	mv *.so ../

clean:
	rm *.so
	touch *.f90
touch-dummy:
	$(dummy-variable)
