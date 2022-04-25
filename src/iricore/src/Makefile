default:
	f2py -m iri_core -h iri_core.pyf iri_core.f90
	gfortran -c -fPIC iri_core.f90 irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for -O3 -w
	#f2py -c iri_core.pyf iri_core.o irisub.o irifun.o iritec.o iridreg.o igrf.o cira.o iriflip.o -DF2PY_REPORT_ATEXIT
	f2py -c iri_core.pyf iri_core.o irisub.o irifun.o iritec.o iridreg.o igrf.o cira.o iriflip.o
	rm ./*.pyf ./*.o ./*.mod
	mv iri_core*.so ../