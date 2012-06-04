      program zahn_bubbles

c     THIS PROGRAM CALCULATES A CATALOG OF IONIZED BUBBLES
c     USING THE METHOD OF ZAHN ET AL. (2006)
c     THE INPUT IS IONIZED FRACTION ON A MESH, THE OUTPUT 
c     IS THE NUMBER OF CELLS OF EACH SIZE 

      character *128 filein1,filein2,fileout
      parameter (pi=3.14159)

      complex, allocatable:: forig(:,:,:),fc1(:,:,:)
      real, allocatable:: x(:,:,:),pk1(:),sizes(:,:,:,:),xs(:,:,:),
     &xin(:,:,:)
      integer, allocatable:: inbin(:,:)
      double precision average1, average2, var1, var2

      if(iargc().le.1) then
         write(*,99)
 99      format('usage: zahn_bubbles <fin> <fout> ',
     &          '[<n> <L_box> <h> <xth> <nscales>]')
         goto 900
      endif
           
      call getarg(1,filein1)
      call getarg(2,fileout)
      n=i4arg(3,256)
      boxsize=r4arg(4,1.e2)
      h=r4arg(5,7.e-1)
      xth=r4arg(6,9.e-1)
      nscales=i4arg(7,20)

      boxsize=boxsize/h

      n12=n/2
      n21=n12+1

      allocate (x(n,n,n))
      allocate (xin(n,n,n))
      allocate (xs(n,n,n))
      allocate (fc1(n21,n,n))
      allocate (forig(n21,n,n))
      allocate (pk1(n))
      allocate (inbin(2,1000))
      allocate (sizes(2,n,n,n))

      sizes=-1.

c     GET INPUT REAL DATA

c      open(unit=1,file=filein1,form='unformatted')
      open(unit=1,file=filein1,form='binary')
      read(1) n,n,n
      read(1) (((xin(i,j,k),k=1,n),j=1,n),i=1,n)
      close(1)

c     GET MEAN AND CONVERT TO ZERO MEAN FIELD

      average1=0.
      do i=1,n
         do j=1,n
            do k=1,n
               average1=average1+xin(i,j,k)
            enddo
         enddo
      enddo	

      average1=average1/real(n)**3
      do i=1,n
         do j=1,n
            do k=1,n               
               x(i,j,k)=xin(i,j,k)/average1-1.
            enddo
         enddo
      enddo

c     PACK COMPLEX ARRAY WITH DATA CUBE

      do i=1,n12
         do j=1,n
            do k=1,n

               i1=2*i-1
               i2=2*i

               v1=x(i1,j,k)
               v2=x(i2,j,k)
               forig(i,j,k)=cmplx(v1,v2)

            enddo
         enddo
      enddo

c     FFT

      forig(n21,:,:)=0.
      call fft3rks(forig,n,n,n)
      forig=forig/real(n)**3 

c     FILTER THE IONIZED FRACTION ON DIFFERENT SCALES R

      rinit=boxsize
      rfinal=boxsize/real(n)/10.
      dscalel=(log10(rfinal)-log10(rinit))/real(nscales-1)
      dlnscale=-dscalel/log10(exp(1.))

c     LOOP OVER DIFFERENT FILTER SCALES R

      dk=2.*pi/boxsize
      d3k=dk**3
      akmax=dk*n
      pk1=0.
      inbin=0
      itot=0	
      itot1=0
      itot2=0
      Rav1=0.
      Rav2=0.
c      write(*,*) 'dlnscale=', dlnscale
c      write(*,*) 'mine    =', 1.0/(-(nscales-1)/log(1/(10.0*n)))


      open(unit=1,file=fileout)
      open(unit=2,file='center.dat')
      do iscale=1,nscales

         R=boxsize* (1.0/(10.0*n))**((iscale-1.0)/(nscales-1.0))
c         write(*,*)'R1=',R
         R=10.**(log10(rinit)+(iscale-1)*dscalel)
c         write(*,*)'R2=',R
         
         fc1=forig
         do k=1,n
            zk=(k-1)*dk
            if(k.gt.n12) zk=zk-akmax
            
            do j=1,n
               yk=(j-1)*dk
               if(j.gt.n12) yk=yk-akmax

               do i=1,n21
                  xk=(i-1)*dk
                  if(i.gt.n12) xk=xk-akmax
                  
                  rk=sqrt(xk**2+yk**2+zk**2)
                  rkR=rk*R
                  if(rk.ne.0.)
c     &                 fc1(i,j,k)=(W(rkR)**2*forig(i,j,k)) ! FILTER ON SCALE R
     &                 fc1(i,j,k)=(W(rkR)*forig(i,j,k)) ! FILTER ON SCA LE R

               enddo
            enddo
         enddo

c        FFT FROM k TO REAL SPACE

         call fft3krs(fc1,n,n,n)

         average2=0.
         do i=1,n12
            do j=1,n
               do k=1,n
                  i1=2*i-1
                  i2=2*i
                  xs(i1,j,k)=real(fc1(i,j,k))
                  xs(i2,j,k)=aimag(fc1(i,j,k))
                  average2=average2+xs(i1,j,k)+xs(i2,j,k)
               enddo
            enddo
         enddo
c         write(101+iscale,*) n,n
         do i=1,n
            do j=1,n
c               write(101+iscale,*) (xs(i,j,100)+1.)*average1
               do k=1,n
                  xcur=(xs(i,j,k)+1.)*average1
		  xHI=1.-xcur
                  if(xcur.gt.xth.and.sizes(1,i,j,k).lt.0.) then
                     sizes(1,i,j,k)=R
                     inbin(1,iscale)=inbin(1,iscale)+1
		     Rav1=Rav1+R
                  endif
                  if(xHI.gt.xth.and.sizes(2,i,j,k).lt.0.) then
                     sizes(2,i,j,k)=R
                     inbin(2,iscale)=inbin(2,iscale)+1
		     Rav2=Rav2+R
                  endif
               enddo
            enddo
         enddo
         itot1=itot1+inbin(1,iscale)
         itot2=itot2+inbin(2,iscale)

c         write(*,*) 'R, iscale,itot1 = ',R,iscale,itot1,inbin(1,iscale)
c         write(*,*) 'xc = ',(xs(n/2,n/2,n/2)+1)*average1
         write(2,*) R,(xs(n/2,n/2,n/2)+1.)*average1
      enddo

      do iscale=1,nscales
         R=10.**(log10(rinit)+(iscale-1)*dscalel)
         write(1,10) R,real(inbin(1,iscale))/dlnscale/real(itot1),
     &   real(inbin(2,iscale))/dlnscale/real(itot2),inbin(1,iscale),
     &        iscale
 10      format(3(1pe13.6,1x),i8,1x,i5)
      enddo
      close(1)

      write(*,11) filein1,average1,Rav1/real(itot1),Rav2/real(itot2)
 11   format(a128,3(1pe13.6,1x))

      open(unit=1,file='sizes.asci')
      open(unit=2,file='rhoHI.asci')
      write(1,*) n,n
      write(2,*) n,n
      do i=1,n
         do j=1,n
            write(1,*) sizes(1,i,j,n/2)
            write(2,*) xin(i,j,n/2)
         enddo
      enddo

      close(1)
      close(2)

 900  continue

      stop
      end

c==============================================================================

      function W(x)

      W=3.*(sin(x)-x*cos(x))/x**3
c      W=1.0

      return
      end
