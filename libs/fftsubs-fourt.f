cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine FFT13(f,n12,n2,n3,isign)
c  Compute n2*n3 1-D FFTs in place on the first dimension of a complex
c  array of size (n12+1,n2,n3).  Real and imaginary parts are assumed to be
c  in successive storage locations.
	implicit none
	integer, intent(in):: n12,n2,n3,isign
	complex, intent(inout) :: f(n12+1,n2,n3)

	complex, dimension(n12):: temp,work
	integer :: i2, i3
c

c  Do FFTs on 1st dimension.
	do i3 = 1,n3
	  do i2 = 1,n2
	    temp(1:n12)=f(1:n12,i2,i3)
	    call fourt(temp,n12,1,isign,1,work)
	    f(1:n12,i2,i3)=temp(1:n12)
          end do
        end do

	return
	end subroutine FFT13
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine FFT23(f,n1,n2,n3,isign)
c  Compute n1*n3 1-D FFTs in place on the second dimension of a complex
c  array of size (n1,n2,n3).  Real and imaginary parts are assumed to be
c  in successive storage locations.
	implicit none
	integer, intent(in):: n1,n2,n3,isign
	complex, intent(inout) :: f(n1,n2,n3)

	complex, dimension(n2):: temp,work
	integer :: i1, i3
c

c  Do FFTs on 2nd dimension.
	do i3 = 1,n3
	  do i1 = 1,n1
	    temp(1:n2)=f(i1,1:n2,i3)
	    call fourt(temp,n2,1,isign,1,work)
	    f(i1,1:n2,i3)=temp(1:n2)
          end do
        end do

	return
	end subroutine FFT23
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine FFT33(f,n1,n2,n3,isign)
c  Compute n2*n3 1-D FFTs in place on the third dimension of a complex
c  array of size (n1,n2,n3).  Real and imaginary parts are assumed to be
c  in successive storage locations.
	implicit none
	integer, intent(in):: n1,n2,n3,isign
	complex, intent(inout) :: f(n1,n2,n3)

	complex, dimension(n3):: temp,work
	integer :: i1, i2
c

c  Do FFTs on 3rd dimension.
	do i2 = 1,n2
	  do i1 = 1,n1
	    temp(1:n3)=f(i1,i2,1:n3)
	    call fourt(temp,n3,1,isign,1,work)
	    f(i1,i2,1:n3)=temp(1:n3)
          end do
        end do

	return
	end subroutine FFT33
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine FFT12(f,n12,n2,isign)
c  Compute n2 1-D FFTs in place on the first dimension of a complex
c  array of size (n12+1,n2).  Real and imaginary parts are assumed to be
c  in successive storage locations.
	implicit none
	integer, intent(in):: n12,n2,isign
	complex, intent(inout) :: f(n12+1,n2)

	complex, dimension(n12):: temp,work
	integer :: i2
c

c  Do FFTs on 1st dimension.
	do i2 = 1,n2
	  temp(1:n12)=f(1:n12,i2)
	  call fourt(temp,n12,1,isign,1,work)
	  f(1:n12,i2)=temp(1:n12)
        end do

	return
	end subroutine FFT12
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine FFT22(f,n1,n2,isign)
c  Compute n1 1-D FFTs in place on the second dimension of a complex
c  array of size (n1,n2).  Real and imaginary parts are assumed to be
c  in successive storage locations.
	implicit none
	integer, intent(in):: n1,n2,isign
	complex, intent(inout) :: f(n1,n2)

	complex, dimension(n2):: temp,work
	integer :: i1
c

c  Do FFTs on 2nd dimension.
	do i1 = 1,n1
	  temp(1:n2)=f(i1,1:n2)
	  call fourt(temp,n2,1,isign,1,work)
	  f(i1,1:n2)=temp(1:n2)
        end do

	return
	end subroutine FFT22
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	SUBROUTINE FOURT (DATA,NN,NDIM,ISIGN,IFORM,WORK)
C
C**********************************************************************
C*								      *
C*		THE COOLEY-TUKEY FAST FOURIER TRANSFORM		      *
C*								      *
C*   TRANSFORM(K1,K2,...) = SUM(DATA(J1,J2,...)*EXP(ISIGN*2*SQRT(-1)  *
C*   *((J1-1)*(K1-1)/NN(1)+(J2-1)*(K2-1)/NN(2)+...))),SUMMED FOR ALL  *
C*   J1, K1 BETWEEN 1 AND NN(1), J2, K2 BETWEEN 1 AND NN(2), ETC.     *
C*   THERE IS NO LIMIT TO THE NUMBER OF SUBSCRIPTS.  DATA IS A        *
C*   MULTIDIMENSIONAL COMPLEX ARRAY WHOSE REAL AND IMAGINARY	      *
C*   PARTS ARE ADJACENT IN STORAGE, SUCH AS FORTRAN IV PLACES THEM.   *
C*   IF ALL IMAGINARY PARTS ARE ZERO (DATA ARE DISGUISED REAL), SET   *
C*   IFORM TO ZERO TO CUT THE RUNNING TIME BY UP TO FORTY PERCENT.    *
C*   OTHERWISE, IFORM = +1. THE LENGTHS OF ALL DIMENSIONS ARE	      *
C*   STORED IN ARRAY NN, OF LENGTH NDIM. THEY MAY BE ANY POSITIVE     *
C*   INTEGERS, THO THE PROGRAM RUNS FASTER ON COMPOSITE INTEGERS, AND *
C*   ESPECIALLY FAST ON NUMBERS RICH IN FACTORS OF TWO. ISIGN IS +1   *
C*   OR -1.  IF A -1 TRANSFORM IS FOLLOWED BY A +1 ONE (OR A +1	      *
C*   BY A -1) THE ORIGINAL DATA REAPPEAR, MULTIPLIED BY NTOT (NN(1)   *
C*   *NN(2)*...).TRANSFORM VALUES ARE ALWAYS COMPLEX,AND ARE RETURNED *
C*   IN ARRAY DATA, REPLACING THE INPUT. IN ADDITION, IF ALL	      *
C*   DIMENSIONS ARE NOT POWERS OF TWO, ARRAY WORK MUST BE SUPPLIED,   *
C*   COMPLEX OF LENGTH EQUAL TO THE LARGEST NON 2**K DIMENSION.       *
C*   OTHERWISE, REPLACE WORK BY ZERO IN THE CALLING SEQUENCE.	      *
C*   NORMAL FORTRAN DATA ORDERING IS EXPECTED, FIRST SUBSCRIPT	      *
C*   VARYING FASTEST. ALL SUBSCRIPTS BEGIN AT ONE.		      *
C*								      *
C*   RUNNING TIME IS MUCH SHORTER THAN THE NAIVE NTOT**2, BEING       *
C*   GIVEN BY THE FOLLOWING FORMULA.  DECOMPOSE NTOT INTO	      *
C*   2**K2 * 3**K3 * 5**K5 * ....LET SUM2=2*K2,SUMF=3*K3+5*K5+...     *
C*   AND NF=K3+K5+...  THE TIME TAKEN BY A MULTIDIMENSIONAL	      *
C*   TRANSFORM ON THESE NTOT DATA IS T = T0 + NTOT*(T1+T2*SUM2+       *
C*   +T3*SUMF+T4*NF). ON THE CDC 3300 (FLOATING POINT ADD TIME OF     *
C*   SIX MICROSECONDS), T = 3000 + NTOT * (500+43*SUM2+68*SUMF+       *
C*   +320*NF) MICROSECONDS ON COMPLEX DATA. IN ADDITION, THE	      *
C*   ACCURACY IS GREATLY IMPROVED, AS THE RMS RELATIVE ERROR IS       *
C*   BOUNDED BY 3*2**(-B)*SUM(FACTOR(J)**1.5), WHERE B IS THE NUMBER  *
C*   OF BITS IN THE FLOATING POINT FRACTION AND FACTOR(J) ARE THE     *
C*   PRIME FACTORS OF NTOT.					      *
C*								      *
C*   THE DISCRETE FOURIER TRANSFORM PLACES THREE RESTRICTIONS UPON    *
C*   THE DATA:							      *
C*   1.  THE NUMBER OF INPUT DATA AND THE NUMBER OF TRANSFORM VALUES  *
C*	 MUST BE THE SAME.					      *
C*   2.  BOTH THE INPUT DATA AND THE TRANSFORM VALUES MUST REPRESENT  *
C*	 EQUISPACED POINTS IN THEIR RESPECTIVE DOMAINS OF TIME AND    *
C*	 FREQUENCY.CALLING THESE SPACINGS DELTAT AND DELTAF, IT MUST  *
C*	 BE TRUE THAT DELTAF=2*PI/(NN(I)*DELTAT).OF COURSE, DELTAT    *
C*	 NEED NOT BE THE SAME FOR EVERY DIMENSION.		      *
C*   3.  CONCEPTUALLY AT LEAST, THE INPUT DATA AND THE TRANSFORM      *
C*	 OUTPUT REPRESENT SINGLE CYCLES OF PERIODIC FUNCTIONS.	      *
C*								      *
C*	DC-COMPONENT IS MULTIPLIED BY NN, OTHERS BY NN/2	      *
C*								      *
C*   EXAMPLE 1. THREE-DIMENSIONAL FORWARD FFT OF A COMPLEX ARRAY      *
C*		DIMENSIONED 32 BY 25 BY 13 IN FORTRAN IV.	      *
C*		DIMENSION DATA(32,25,13),WORK(50),NN(3)		      *
C*		COMPLEX   DATA					      *
C*		DATA	  NN	/32,25,13/			      *
C*		DO 1 I = 1,32					      *
C*		DO 1 J = 1,25					      *
C*		DO 1 K = 1,13					      *
C*   1		DATA(I,J,K) = COMPLEX VALUE			      *
C*		CALL FOURT (DATA,NN,3,-1,1,WORK)		      *
C*								      *
C*   EXAMPLE 2. ONE-DIMENSIONAL FORWARD FFT OF A REAL ARRAY OF	      *
C*		LENGTH 64 IN FORTRAN II.			      *
C*		DIMENSION DATA(2,64)				      *
C*		DO 2 I = 1,64					      *
C*		DATA(1,I) = REAL PART				      *
C*   2		DATA(2,I) = 0.					      *
C*		CALL FOURT (DATA,64,1,-1,0,0)			      *
C*								      *
C**********************************************************************
C
	DIMENSION DATA(*),NN(*),IFACT(32),WORK(*)
	TWOPI = 6.2831853
	IF (NDIM-1) 920,1,1
    1	NTOT = 2
	DO 2 IDIM = 1,NDIM
	  IF (NN(IDIM)) 920,920,2
    2	NTOT = NTOT * NN(IDIM)
C
C		*** MAIN LOOP FOR EACH DIMENSION ***
C
	NP1 = 2
	DO 910 IDIM = 1,NDIM
	  N = NN(IDIM)
	  NP2 = NP1 * N
	  IF (N-1) 920,900,5
C
C		*** FACTOR N ***
C
    5	  M = N
	  NTWO = NP1
	  IF = 1
	  IDIV = 2
   10	  IQUOT = M / IDIV
	  IREM = M - IDIV * IQUOT
	  IF (IQUOT-IDIV) 50,11,11
   11	  IF (IREM) 20,12,20
   12	  NTWO = NTWO + NTWO
	  M = IQUOT
	  GO TO 10
   20	  IDIV = 3
   30	  IQUOT = M / IDIV
	  IREM = M - IDIV * IQUOT
	  IF (IQUOT-IDIV) 60,31,31
   31	  IF (IREM) 40,32,40
   32	  IFACT(IF) = IDIV
	  IF = IF + 1
	  M = IQUOT
	  GO TO 30
   40	  IDIV = IDIV + 2
	  GO TO 30
   50	  IF (IREM) 60,51,60
   51	  NTWO = NTWO + NTWO
	  GO TO 70
   60	  IFACT(IF) = M
C
C	SEPARATE FOUR CASES--
C	  1. COMPLEX TRANSFORM OR REAL TRANSFORM FOR THE 4TH,5TH,ETC.
C	     DIMENSIONS.
C	  2. REAL TRANSFORM FOR THE 2ND OR 3RD DIMENSION.  METHOD--
C	     TRANSFORM HALF THE DATA, SUPPLYING THE OTHER HALF BY
C	     CONJUGATE SYMMETRY.
C	  3. REAL TRANSFORM FOR THE 1ST DIMENSION, N ODD.  METHOD--
C	     TRANSFORM HALF THE DATA AT EACH STAGE, SUPPLYING THE
C	     OTHER HALF BY CONJUGATE SYMMETRY.
C	  4. REAL TRANSFORM FOR THE 1ST DIMENSION,N EVEN.  METHOD--
C	     TRANSFORM A COMPLEX ARRAY OF LENGTH N/2 WHOSE REAL PARTS
C	     ARE THE EVEN-NUMBERED REAL VALUES AND WHOSE IMAGINARY
C	     PARTS ARE THE ODD-NUMBERED REAL VALUES. SEPARATE AND
C	     SUPPLY THE SECOND HALF BY CONJUGATE SYMMETRY.
C
   70	  NON2 = NP1 * (NP2 / NTWO)
	  ICASE = 1
	  IF (IDIM-4) 71,90,90
   71	  IF (IFORM) 72,72,90
   72	  ICASE = 2
	  IF (IDIM-1) 73,73,90
   73	  ICASE = 3
	  IF (NTWO-NP1) 90,90,74
   74	  ICASE = 4
	  NTWO = NTWO / 2
	  N = N / 2
	  NP2 = NP2 / 2
	  NTOT = NTOT / 2
	  I = 3
	  DO 80 J = 2,NTOT
	    DATA(J) = DATA(I)
	    I = I + 2
   80	  CONTINUE
   90	  I1RNG = NP1
	  IF (ICASE-2) 100,95,100
   95	  I1RNG = NP0 * (1 + NPREV / 2)
C
C	*** SHUFFLE ON THE FACTORS OF TWO IN N. AS THE SHUFFLING
C	*** CAN BE DONE BY SIMPLE INTERCHANGE, NO WORK AREA IS NEEDED.
C
  100	  IF (NTWO-NP1) 600,600,110
  110	  NP2HF = NP2 / 2
	  J = 1
	  DO 150 I2 = 1,NP2,NON2
	    IF (J-I2) 120,130,130
  120	    I1MAX = I2 + NON2 - 2
	    DO 125 I1 = I2,I1MAX,2
	      DO 125 I3 = I1,NTOT,NP2
	        J3 = J + I3 - I2
	        TEMPR = DATA(I3)
		TEMPI = DATA(I3+1)
		DATA(I3) = DATA(J3)
		DATA(I3+1) = DATA(J3+1)
		DATA(J3) = TEMPR
		DATA(J3+1) = TEMPI
  125	    CONTINUE
  130	    M = NP2HF
  140	    IF (J-M) 150,150,145
  145	    J = J - M
	    M = M / 2
	    IF (M-NON2) 150,140,140
  150	  J = J + M
C
C	MAIN LOOP FOR FACTORS OF TWO. PERFORM FOURIER TRANSFORMS OF
C	LENGTH FOUR, WITH ONE OF LENGTH TWO IF NEEDED. THE
C	TWIDDLE FACTOR W=EXP(ISIGN*2*PI*SQRT(-1)*M/(4*MMAX)). CHECK
C	FOR W=ISIGN*SQRT(-1) AND REPEAT FOR W=SQRT(-1)*CONJUGATE(W).
C
	  NON2T = NON2 + NON2
	  IPAR = NTWO / NP1
  310	  IF (IPAR-2) 350,330,320
  320	  IPAR = IPAR / 4
	  GO TO 310
  330	  CONTINUE
	  DO 340 I1 = 1,I1RNG,2
	    DO 340 J3 = I1,NON2,NP1
	      DO 340 K1 = J3,NTOT,NON2T
		K2 = K1 + NON2
		TEMPR = DATA(K2)
		TEMPI = DATA(K2+1)
		DATA(K2) = DATA(K1) - TEMPR
		DATA(K2+1) = DATA(K1+1) - TEMPI
		DATA(K1) = DATA(K1) + TEMPR
		DATA(K1+1) = DATA(K1+1) + TEMPI
  340	  CONTINUE
  350	  MMAX = NON2
  360	  IF (MMAX-NP2HF) 370,600,600
  370	  LMAX = MAX0 (NON2T,MMAX/2)
	  IF (MMAX-NON2) 405,405,380
  380	  THETA = -TWOPI * FLOAT(NON2)/FLOAT(4*MMAX)
	  IF (ISIGN) 400,390,390
  390	  THETA = -THETA
  400	  WR = COS (THETA)
	  WI = SIN (THETA)
	  WSTPR = -2. * WI * WI
	  WSTPI =  2. * WR * WI
  405	  CONTINUE
	  DO 570 L = NON2,LMAX,NON2T
	    M = L
	    IF (MMAX-NON2) 420,420,410
  410	    W2R = WR * WR - WI * WI
	    W2I = 2. * WR * WI
	    W3R = W2R * WR - W2I * WI
	    W3I = W2R * WI + W2I * WR
  420	   CONTINUE
	    DO 530 I1 = 1,I1RNG,2
	      DO 530 J3 = I1,NON2,NP1
		KMIN = J3 + IPAR * M
		IF (MMAX-NON2) 430,430,440
  430		KMIN = J3
  440		KDIF = IPAR * MMAX
  450		KSTEP = 4 * KDIF
		DO 520 K1 = KMIN,NTOT,KSTEP
		  K2 = K1 + KDIF
		  K3 = K2 + KDIF
		  K4 = K3 + KDIF
		  IF (MMAX-NON2) 460,460,480
  460		  U1R = DATA(K1) + DATA(K2)
		  U1I = DATA(K1+1) + DATA(K2+1)
		  U2R = DATA(K3) + DATA(K4)
		  U2I = DATA(K3+1) + DATA(K4+1)
		  U3R = DATA(K1) - DATA(K2)
		  U3I = DATA(K1+1) - DATA(K2+1)
		  IF (ISIGN) 470,475,475
  470		  U4R = DATA(K3+1) -DATA(K4+1)
		  U4I = DATA(K4) - DATA(K3)
		  GO TO 510
  475		  U4R = DATA(K4+1) - DATA(K3+1)
		  U4I = DATA(K3) - DATA(K4)
		  GO TO 510
  480		  T2R = W2R * DATA(K2) - W2I * DATA(K2+1)
		  T2I = W2R * DATA(K2+1) + W2I * DATA(K2)
		  T3R = WR * DATA(K3) - WI * DATA(K3+1)
		  T3I = WR * DATA(K3+1) + WI * DATA(K3)
		  T4R = W3R * DATA(K4) - W3I * DATA(K4+1)
		  T4I = W3R * DATA(K4+1) + W3I * DATA(K4)
		  U1R = DATA(K1) + T2R
		  U1I = DATA(K1+1) + T2I
		  U2R = T3R + T4R
		  U2I = T3I + T4I
		  U3R = DATA(K1) - T2R
		  U3I = DATA(K1+1) - T2I
		  IF (ISIGN) 490,500,500
  490		  U4R = T3I - T4I
		  U4I = T4R - T3R
		  GO TO 510
  500		  U4R = T4I - T3I
		  U4I = T3R - T4R
  510		  DATA(K1) = U1R + U2R
		  DATA(K1+1) = U1I + U2I
		  DATA(K2) = U3R + U4R
		  DATA(K2+1) = U3I + U4I
		  DATA(K3) = U1R - U2R
		  DATA(K3+1) = U1I - U2I
		  DATA(K4) = U3R - U4R
		  DATA(K4+1) = U3I - U4I
  520		CONTINUE
		KMIN = 4 * (KMIN-J3) + J3
		KDIF = KSTEP
		IF (KDIF-NP2) 450,530,530
  530	      CONTINUE
	      M = MMAX - M
	      IF (ISIGN) 540,550,550
  540	      TEMPR = WR
	      WR = -WI
	      WI = -TEMPR
	      GO TO 560
  550	      TEMPR = WR
	      WR = WI
	      WI = TEMPR
  560	      IF (M-LMAX) 565,565,410
  565	      TEMPR = WR
	      WR = WR * WSTPR - WI * WSTPI + WR
	      WI = WI * WSTPR + TEMPR * WSTPI + WI
  570	  CONTINUE
	  IPAR = 3 - IPAR
	  MMAX = MMAX + MMAX
	  GO TO 360
C
C	MAIN LOOP FOR FACTORS NOT EQUAL TO TWO. APPLY THE TWIDDLE
C	FACTOR W=EXP(ISIGN*2*PI*SQRT(-1)*(J2-1)*(J1-J2)/(NP2*IFP1)),
C	THEN PERFORM A FOURIER TRANSFORM OF LENGTH IFACT(IF),MAKING
C	USE OF CONJUGATE SYMMETRIES.
C
  600	  IF (NTWO-NP2) 605,700,700
  605	  IFP1 = NON2
	  IF = 1
	  NP1HF = NP1 / 2
  610	  IFP2 = IFP1 / IFACT(IF)
	  J1RNG = NP2
	  IF (ICASE-3) 612,611,612
  611	  J1RNG = (NP2+IFP1) / 2
	  J2STP = NP2 / IFACT(IF)
	  J1RG2 = (J2STP+IFP2) / 2
  612	  J2MIN = 1 + IFP2
	  IF (IFP1-NP2) 615,640,640
  615	  CONTINUE
	  DO 635 J2 = J2MIN,IFP1,IFP2
	    THETA = -TWOPI * FLOAT(J2-1) / FLOAT(NP2)
	    IF (ISIGN) 625,620,620
  620	    THETA = -THETA
  625	    SINTH = SIN (THETA/2.)
	    WSTPR = -2. * SINTH * SINTH
	    WSTPI = SIN (THETA)
	    WR = WSTPR + 1.
	    WI = WSTPI
	    J1MIN = J2 + IFP1
	    DO 635 J1 = J1MIN,J1RNG,IFP1
	      I1MAX = J1 + I1RNG - 2
	      DO 630 I1 = J1,I1MAX,2
		DO 630 I3 = I1,NTOT,NP2
		  J3MAX = I3 + IFP2 - NP1
		  DO 630 J3 = I3,J3MAX,NP1
		    TEMPR = DATA(J3)
		    DATA(J3) = DATA(J3) * WR - DATA(J3+1) * WI
		    DATA(J3+1) = TEMPR * WI + DATA(J3+1) * WR
  630	      CONTINUE
	      TEMPR = WR
	      WR = WR * WSTPR - WI * WSTPI + WR
	      WI = TEMPR * WSTPI + WI * WSTPR + WI
  635	  CONTINUE
  640	  THETA = -TWOPI / FLOAT(IFACT(IF))
	  IF (ISIGN) 650,645,645
  645	  THETA = -THETA
  650	  SINTH = SIN (THETA/2.)
	  WSTPR = -2. * SINTH * SINTH
	  WSTPI = SIN (THETA)
	  KSTEP = 2 * N / IFACT(IF)
	  KRANG = KSTEP * (IFACT(IF)/2) + 1
	  DO 698 I1 = 1,I1RNG,2
	    DO 698 I3 = I1,NTOT,NP2
	      DO 690 KMIN = 1,KRANG,KSTEP
		J1MAX = I3 + J1RNG - IFP1
		DO 680 J1 = I3,J1MAX,IFP1
		  J3MAX = J1 + IFP2 - NP1
		  DO 680 J3 = J1,J3MAX,NP1
		    J2MAX = J3 + IFP1 - IFP2
		    K = KMIN + (J3-J1+(J1-I3)/IFACT(IF)) / NP1HF
		    IF (KMIN-1) 655,655,665
  655		    SUMR = 0.
		    SUMI = 0.
		    DO 660 J2 = J3,J2MAX,IFP2
		      SUMR = SUMR + DATA(J2)
		      SUMI = SUMI + DATA(J2+1)
  660		    CONTINUE
		    WORK(K) = SUMR
		    WORK(K+1) = SUMI
		    GO TO 680
  665		    KCONJ = K + 2 * (N-KMIN+1)
		    J2 = J2MAX
		    SUMR = DATA(J2)
		    SUMI = DATA(J2+1)
		    OLDSR = 0.
		    OLDSI = 0.
		    J2 = J2 - IFP2
  670		    TEMPR = SUMR
		    TEMPI = SUMI
		    SUMR = TWOWR * SUMR - OLDSR + DATA(J2)
		    SUMI = TWOWR * SUMI - OLDSI + DATA(J2+1)
		    OLDSR = TEMPR
		    OLDSI = TEMPI
		    J2 = J2 - IFP2
		    IF (J2-J3) 675,675,670
  675		    TEMPR = WR * SUMR - OLDSR + DATA(J2)
		    TEMPI = WI * SUMI
		    WORK(K) = TEMPR - TEMPI
		    WORK(KCONJ) = TEMPR + TEMPI
		    TEMPR = WR * SUMI - OLDSI + DATA(J2+1)
		    TEMPI = WI * SUMR
		    WORK(K+1) = TEMPR + TEMPI
		    WORK(KCONJ+1) = TEMPR - TEMPI
  680		CONTINUE
		IF (KMIN-1) 685,685,686
  685		WR = WSTPR + 1.
		WI = WSTPI
		GO TO 690
  686		TEMPR = WR
		WR = WR * WSTPR - WI * WSTPI + WR
		WI = TEMPR * WSTPI + WI * WSTPR + WI
  690	      TWOWR = WR + WR
	      IF (ICASE-3) 692,691,692
  691	      IF (IFP1-NP2) 695,692,692
  692	      K = 1
	      I2MAX = I3 + NP2 - NP1
	      DO 693 I2 = I3,I2MAX,NP1
		DATA(I2) = WORK(K)
		DATA(I2+1) = WORK(K+1)
		K = K + 2
  693	      CONTINUE
	      GO TO 698
C
C	COMPLETE A REAL TRANSFORM IN THE 1ST DIMENSION, N ODD, BY
C	CONJUGATE SYMMETRIES AT EACH STAGE.
C
  695	      J3MAX = I3 + IFP2 - NP1
	      DO 697 J3 = I3,J3MAX,NP1
		J2MAX = J3 + NP2 - J2STP
		DO 697 J2 = J3,J2MAX,J2STP
		  J1MAX = J2 + J1RG2 - IFP2
		  J1CNJ = J3 + J2MAX + J2STP - J2
		  DO 697 J1 = J2,J1MAX,IFP2
		    K = 1 + J1 - I3
		    DATA(J1) = WORK(K)
		    DATA(J1+1) = WORK(K+1)
		    IF (J1-J2) 697,697,696
  696		    DATA(J1CNJ) = WORK(K)
		    DATA(J1CNJ+1) = -WORK(K+1)
  697	      J1CNJ = J1CNJ - IFP2
  698	    CONTINUE
	    IF = IF + 1
	    IFP1 = IFP2
	    IF (IFP1-NP1) 700,700,610
C
C	COMPLETE A REAL TRANSFORM IN THE 1ST DIMENSION, N EVEN, BY
C	CONJUGATE SYMMETRIES.
C
  700	  GO TO (900,800,900,701) , ICASE
  701	  NHALF = N
	  N = N + N
	  THETA = -TWOPI / FLOAT(N)
	  IF (ISIGN) 703,702,702
  702	  THETA = -THETA
  703	  SINTH = SIN (THETA/2.)
	  WSTPR = -2. * SINTH * SINTH
	  WSTPI = SIN (THETA)
	  WR = WSTPR + 1.
	  WI = WSTPI
	  IMIN = 3
	  JMIN = 2 * NHALF - 1
	  GO TO 725
  710	  J = JMIN
	  DO 720 I = IMIN,NTOT,NP2
	    SUMR = (DATA(I)+DATA(J)) / 2.
	    SUMI = (DATA(I+1)+DATA(J+1)) / 2.
	    DIFR = (DATA(I)-DATA(J)) / 2.
	    DIFI = (DATA(I+1)-DATA(J+1)) / 2.
	    TEMPR = WR * SUMI + WI * DIFR
	    TEMPI = WI * SUMI - WR * DIFR
	    DATA(I) = SUMR + TEMPR
	    DATA(I+1) = DIFI + TEMPI
	    DATA(J) = SUMR - TEMPR
	    DATA(J+1) = -DIFI + TEMPI
  720	  J = J + NP2
	  IMIN = IMIN + 2
	  JMIN = JMIN - 2
	  TEMPR = WR
	  WR = WR * WSTPR - WI * WSTPI + WR
	  WI = TEMPR * WSTPI + WI * WSTPR + WI
  725	  IF (IMIN-JMIN) 710,730,740
  730	  IF (ISIGN) 731,740,740
  731	  CONTINUE
	  DO 735 I = IMIN,NTOT,NP2
  735	  DATA(I+1) = -DATA(I+1)
  740	  NP2 = NP2 + NP2
	  NTOT = NTOT + NTOT
	  J = NTOT + 1
	  IMAX = NTOT / 2 + 1
  745	  IMIN = IMAX - 2 * NHALF
	  I = IMIN
	  GO TO 755
  750	  DATA(J) = DATA(I)
	  DATA(J+1) = -DATA(I+1)
  755	  I = I + 2
	  J = J - 2
	  IF (I-IMAX) 750,760,760
  760	  DATA(J) = DATA(IMIN) - DATA(IMIN+1)
	  DATA(J+1) = 0.
	  IF (I-J) 770,780,780
  765	  DATA(J) = DATA(I)
	  DATA(J+1) = DATA(I+1)
  770	  I = I - 2
	  J = J - 2
	  IF (I-IMIN) 775,775,765
  775	  DATA(J) = DATA(IMIN) + DATA(IMIN+1)
	  DATA(J+1) = 0.
	  IMAX = IMIN
	  GO TO 745
  780	  DATA(1) = DATA(1) + DATA(2)
	  DATA(2) = 0.
	  GO TO 900
C
C	COMPLETE A REAL TRANSFORM FOR THE 2ND OR 3RD DIMENSION
C	BY CONJUGATE SYMMETRIES.
C
  800	  IF (I1RNG-NP1) 805,900,900
  805	  CONTINUE
	  DO 860 I3 = 1,NTOT,NP2
	    I2MAX = I3 + NP2 - NP1
	    DO 860 I2 = I3,I2MAX,NP1
	      IMIN = I2 + I1RNG
	      IMAX = I2 + NP1 - 2
	      JMAX = 2 * I3 + NP1 - IMIN
	      IF (I2-I3) 820,820,810
  810	      JMAX = JMAX + NP2
  820	      IF (IDIM-2) 850,850,830
  830	      J = JMAX + NP0
	      DO 840 I = IMIN,IMAX,2
		DATA(I) = DATA(J)
		DATA(I+1) = -DATA(J+1)
  840	      J = J - 2
  850	      J = JMAX
	      DO 860 I = IMIN,IMAX,NP0
		DATA(I) = DATA(J)
		DATA(I+1) = -DATA(J+1)
  860	  J = J - NP0
C
C	END OF LOOP FOR EACH DIMENSION
C
  900	  NP0 = NP1
	  NP1 = NP2
  910	NPREV = N
  920	RETURN
	END
