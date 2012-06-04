F77 	    = ifort
F77_OPT     = -O3 
F77_COMPILE = $(F77_OPT) -c
LIBDIR      = /home/martina/lib/ANALYSIS_TOOLS

.f.o: 
	$(F77) $(F77_COMPILE) $*.f 

ffttest: ffttest.o
	$(F77) $(F77_OPT) -L$(LIBDIR) ffttest.o -lfft -ltext -o ffttest
	cp ffttest ~/bin
pk_spherical: pk_spherical.o
	$(F77) $(F77_OPT) -L$(LIBDIR) pk_spherical.o -lfft -ltext -o pk_spherical
	cp pk_spherical ~/bin
cc_spherical: cc_spherical.o
	$(F77) $(F77_OPT) -L$(LIBDIR) cc_spherical.o -lfft -ltext -o cc_spherical
	cp cc_spherical ~/bin
cc_all: cc_all.o
	$(F77) $(F77_OPT) -L$(LIBDIR) cc_all.o -lfft -ltext -o cc_all
	cp cc_all ~/bin
clean:
	rm *.o 
