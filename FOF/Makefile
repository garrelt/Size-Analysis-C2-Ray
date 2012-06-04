F77 	    = ifort
F77_OPT     = -O3 
F77_COMPILE = $(F77_OPT) -c
LIBDIR      = /home/martina/lib/ANALYSIS_TOOLS

.f.o: 
	$(F77) $(F77_COMPILE) $*.f 

fof_bubbles:  fof_bubbles.o
	$(F77) $(F77_OPT) -L$(LIBDIR) fof_bubbles.o -ltext -o fof_bubbles
	cp fof_bubbles $(HOME)/bin/
fof_bubbles64:  fof_bubbles.o
	$(F77) $(F77_OPT) -L$(LIBDIR) fof_bubbles.o -ltext64 -o fof_bubbles64
	cp fof_bubbles64 $(HOME)/bin/
dndm_bubbles:  dndm_bubbles.o
	$(F77) $(F77_OPT) -L$(LIBDIR) dndm_bubbles.o -ltext -o dndm_bubbles
	cp dndm_bubbles $(HOME)/bin/

fof_bubbles_correlations: fof_bubbles_correlations.o
	$(F77) $(F77_OPT) -L$(LIBDIR) fof_bubbles_correlations.o -ltext -o fof_bubbles_correlations
	cp fof_bubbles_correlations $(HOME)/bin/

clean:
	rm *.o *.f~ *.inc~ *.anl *.m *.list
