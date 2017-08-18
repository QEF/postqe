
QE_DIR = /home/brunato/Downloads/QE_Forge/QE/q-e
#QE_DIR = /home/mauropalumbo/q-e

.SUFFIXES :
.SUFFIXES : .o .c .f .f90


.f90.o:
	gfortran -c -fPIC $<

pyqe:   wrappers_utils.f90 wrappers_funct.f90 
	f2py3 -c wrappers_utils.f90 wrappers_funct.f90 -m pyqe -I$(QE_DIR)/Modules/ -I$(QE_DIR)/clib/ \
	-llapack -lblas libqemod.a clib.a libqefft.a libpw.a

pyqe2:  wrappers_vloc.f90
	f2py3 -c wrappers_vloc.f90 -m pyqe2 libQE.a
	
all: pyqe pyqe2 move

move :
	mv *.so ../

clean:
	rm *.so
	touch *.f90
touch-dummy:
	$(dummy-variable)
