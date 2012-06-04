# standard flags for AIX
CC	= gcc -w
LD	= gcc
COPT	= -O3
CFLAGS	= $(COPT)
LFLAGS	= -lm -s
TIMEX	= timex
MAIN	= minkowski

# additional flags for HP-UX on HP 7 series
#CC	= gcc -Aa
#COPT 	= -O4
#LFLAGS	= -lf -lm -s
# additional flags for a Dec alpha
#CC	= gcc
#LD	= gcc
#TIMEX	= time
# additional flags for SGIs
#LFLAGS	= -lm -lmalloc -s
# additional flags for SunOS on Sun Sparc Stations
#CC	= gcc
#LD	= gcc
#TIMEX	= /usr/5bin/time

OBJECTS	= fourier.o image.o main.o minkowski_th.o org.o
#OBJECTS	= fourier.o image.o main.o minkowski.o org.o
$(MAIN) main:	$(OBJECTS)
		$(LD) $(OBJECTS) $(LFLAGS) -o $(MAIN)
		cp $(MAIN) $(HOME)/bin

test:		$(MAIN)
		./$(MAIN) -\?
		$(TIMEX) ./$(MAIN)\
		-l-4 -h4 -b100 -s4 -L10 -x32 -y32 -z32 -otest.res -0test.xpm
# xv -geometry +0+0 test.xpm &
.c.o:		
		$(CC) -c $(CFLAGS) $<
.f.o:		
		$(FC) -c $(FFLAGS) $<

clean:	
		-rm -f *.o *~ \#* core

#realclean:	clean
#		-rm -f test* $(MAIN)

tar:		
		tar cf - *.[cfh] Makefile | gzip --best > source.tar.gz

suicide:	tar realclean
		rm -f *.[cfh] Makefile
