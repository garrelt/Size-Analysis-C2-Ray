F77 	    = ifort
F77_OPT     = -O3 
F77_COMPILE = $(F77_OPT) -c
LIBDIR      = /home/martina/lib/ANALYSIS_TOOLS

.f.o: 
	$(F77) $(F77_COMPILE) $*.f 

zahn_bubbles:  zahn_bubbles.o
	$(F77) $(F77_OPT) -L$(LIBDIR) zahn_bubbles.o -ltext -lfft -o zahn_bubbles
	cp zahn_bubbles $(HOME)/bin/

clean:
	rm *.o *.f~ *.inc~ *.anl *.m *.list
