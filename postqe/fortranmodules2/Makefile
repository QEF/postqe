
.SUFFIXES :
.SUFFIXES : .o .c .f .f90


.f90.o:
	gfortran -c -fPIC $<


libQE.a: touch-dummy
	(cd flib ; make all)

qeutils:  wrappers_utils.f90
	f2py -c wrappers_utils.f90 -m qeutils
	
all: qeutils

move :
	mv *.so ../

clean:
	rm *.so
	touch *.f90
touch-dummy:
	$(dummy-variable)
