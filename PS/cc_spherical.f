      program cc_spherical

      character *128 filein1,filein2,fileout
      parameter (pi=3.14159)

      complex, allocatable:: fc1(:,:,:),fc2(:,:,:), fc3(:,:,:)
      real, allocatable:: data1(:,:,:),data2(:,:,:),data3(:,:,:)
      real, allocatable:: pk1(:),pk2(:), pk3(:),cc(:)
      integer, allocatable:: inbin(:)
      double precision average1, average2, average3, var1, var2

      if(iargc().le.1) then
         write(*,99)
 99      format('usage: cc_spherical <f1> <f2> <fout> ',
     &          '[<n> <Lbox> <h>]')
         goto 900
      endif
           
      call getarg(1,filein1)
      call getarg(2,filein2)
      call getarg(3,fileout)
      n=i4arg(4,203)
      boxsize=r4arg(5,1.e2)
      h=r4arg(6,7.e-1)

      n12=n/2
      n21=n12+1

      allocate (data1(n,n,n))
      allocate (data2(n,n,n))
      allocate (data3(n,n,n))
      allocate (fc1(n21,n,n))
      allocate (fc2(n21,n,n))
      allocate (fc3(n21,n,n))
      allocate (pk1(n))
      allocate (pk2(n))
      allocate (pk3(n))
      allocate (cc(n))
      allocate (inbin(n))

c     GET INPUT REAL DATA

      open(unit=1,file=filein1,form='binary')
      read(1) n,n,n
      read(1) (((data1(i,j,k),k=1,n),j=1,n),i=1,n)
      close(1)
      open(unit=1,file=filein2,form='binary')
      read(1) n,n,n
      write(*,*)'n=',n
      read(1) (((data2(i,j,k),k=1,n),j=1,n),i=1,n)
      close(1)
      data3=data1*(1.0-data2)

c     GET MEAN AND CONVERT TO ZERO MEAN FIELD

      average1=0.
      average2=0.
      average3=0.
      do i=1,n
         do j=1,n
            do k=1,n
               average1=average1+data1(i,j,k)
               average2=average2+data2(i,j,k)
               average3=average3+data3(i,j,k)
            enddo
         enddo
      enddo

      average1=average1/real(n)**3
      average2=average2/real(n)**3
      average3=average3/real(n)**3

      data1=data1/average1-1.
      data2=data2-average2 ! we do not normalize here
      data3=data3/average3-1.

c     PACK COMPLEX ARRAY WITH DATA CUBE

      do i=1,n12
         do j=1,n
            do k=1,n

               i1=2*i-1
               i2=2*i

               v1=data1(i1,j,k)
               v2=data1(i2,j,k)
               fc1(i,j,k)=cmplx(v1,v2)

               v1=data2(i1,j,k)
               v2=data2(i2,j,k)
               fc2(i,j,k)=cmplx(v1,v2)

               v1=data3(i1,j,k)
               v2=data3(i2,j,k)
               fc3(i,j,k)=cmplx(v1,v2)

            enddo
         enddo
      enddo

c     FFT

      fc1(n21,:,:)=0.
      fc2(n21,:,:)=0.
      fc3(n21,:,:)=0.
      call fft3rks(fc1,n,n,n)
      call fft3rks(fc2,n,n,n)
      call fft3rks(fc3,n,n,n)
      fc1=fc1/real(n)**3 
      fc2=fc2/real(n)**3 
      fc3=fc3/real(n)**3 
c     SPHERICAL AVERAGE

      boxsize=boxsize/h
      dk=2.*pi/boxsize
      d3k=dk**3
      akmax=dk*n

      pk1=0.
      pk2=0.
      pk3=0.
      cc=0.
      inbin=0

      var1=0.
      var2=0.

      do k=1,n
         zk=(k-1)*dk
         if(k.gt.n12) zk=zk-akmax

         do j=1,n
            yk=(j-1)*dk
            if(j.gt.n12) yk=yk-akmax

            do i=2,n21
               xk=(i-1)*dk
               if(i.gt.n12) xk=xk-akmax

               var1=var1+2.0d0*fc1(i,j,k)*conjg(fc1(i,j,k)) ! ff* = variance per mode
               var2=var2+2.0d0*fc2(i,j,k)*conjg(fc2(i,j,k)) ! ff* = variance per mode

               rk=sqrt(xk**2+yk**2+zk**2)
               ibin=int(sqrt((xk/dk)**2+(yk/dk)**2+(zk/dk)**2))+1

               if(ibin.ge.1.and.ibin.le.n) then

                  pcur1=fc1(i,j,k)*conjg(fc1(i,j,k))/dk**3*(2.*pi)**3 ! Delta(k)=pcur*k^3/(2pi^2)
                  pcur2=fc2(i,j,k)*conjg(fc2(i,j,k))/dk**3*(2.*pi)**3 ! Delta(k)=pcur*k^3/(2pi^2)
                  pcur3=fc3(i,j,k)*conjg(fc3(i,j,k))/dk**3*(2.*pi)**3 ! Delta(k)=pcur*k^3/(2pi^2)
                  cccur=fc1(i,j,k)*conjg(fc2(i,j,k))/dk**3*(2.*pi)**3 

                  inbin(ibin)=inbin(ibin)+1
                  pk1(ibin)=pk1(ibin)+pcur1
                  pk2(ibin)=pk2(ibin)+pcur2
                  pk3(ibin)=pk3(ibin)+pcur3
                  cc(ibin)=cc(ibin)+cccur

               endif

            enddo
            var1=var1+fc1(1,j,k)*conjg(fc1(1,j,k))
            var2=var2+fc2(1,j,k)*conjg(fc2(1,j,k))
         enddo
      enddo

      open(unit=1,file=fileout)

      var1=0.
      var2=0.
      pkmax=0.
      do i=1,n

         pk1(i)=pk1(i)/real(inbin(i))
         pk2(i)=pk2(i)/real(inbin(i))
         pk3(i)=pk3(i)/real(inbin(i))
         cc(i)=cc(i)/real(inbin(i))

	 fkout=dk*real(i)
	 fkin=dk*real(i-1)
	 fkcur=dk*(real(i)-0.5)
         d3k=4.*pi/3.*(fkout**3-fkin**3)

         if(inbin(i).gt.0) then
            var1=var1+d3k*pk1(i)/(8.*pi**3) ! variance
            var2=var2+d3k*pk2(i)/(8.*pi**3) ! variance
            write(1,10) fkcur,pk1(i),pk2(i),cc(i), pk3(i),inbin(i)
	    dcur=pk2(i)*fkcur**3
	    if(dcur.gt.pkmax) then
		peak_xx=fkcur
		pkmax=dcur
	    endif
         endif
 10      format(5(1pe13.6,1x),i6.6)
      enddo

      close(1)

      write(*,11) filein1,average2,peak_xx
 11   format(a128,1x,2(1pe13.6,1x))

 900  continue

      stop
      end
