
celldm: calcola_celldm.f90
	f2py3 -c -m celldm calcola_celldm.f90
ggen:   generate_gg_list.f90
	f2py3 -c -m ggen generate_gg_list.f90
vloc:	vloc_of_g.f90  simpsn.f90  erf.f90
	f2py3 -c -m vloc vloc_of_g.f90 simpsn.f90 erf.f90
sf:	struct_fact.f90
	f2py3 -c -m sf struct_fact.f90
gshells: gshells.f90
	f2py3 -c -m gshells gshells.f90
funct: funct.f90
	f2py3 -c -m funct funct.f90 functionals.f90 lsda_functionals.f90 expint.f90 erf.f90
all: gshells sf vloc ggen funct celldm

move :
	mv *.so ../

clean:
	rm *.so
	touch *.f90
