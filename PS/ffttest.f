      program fftdriver

      character *128 filien,fileout
      parameter (pi=3.14159)

      complex, allocatable:: fc(:,:,:)
      real, allocatable:: data(:,:,:)
      
      call getarg(1,fileout)
      n=i4arg(2,32)
      f=r4arg(3,5.)

      n12=n/2
      n21=n12+1

      allocate (data(n,n,n))
      allocate (fc(n21,n,n))

      fk=2.*pi*f
      do i=1,n
         x=(real(i)-1.)/real(n-1)
         val=sin(fk*x)
         do j=1,n
            y=(real(j)-0.5)/real(n-1)
            do k=1,n
               z=(real(k)-0.5)/real(n-1)
               data(i,j,k)=val
            enddo
         enddo
      enddo

c     PACK COMPLEX ARRAY WITH DATA CUBE

      do i=1,n12
         do j=1,n
            do k=1,n
               i1=2*(i-1)
               i2=2*(i-1)+1
               v1=data(i1,j,k)
               v2=data(i2,j,k)
               fc(i,j,k)=cmplx(v1,v2)
            enddo
         enddo
      enddo

      fc(n21,:,:)=0.
      call fft3rks(fc,n,n,n)

      dk=2.*pi
      do i=1,n12
         xk=(i-1)*dk
         do j=1,n
            yk=(j-1)*dk
            do k=1,n
               zk=(k-1)*dk
               rk=sqrt(xk**2+yk**2+zk**2)
               pcur=abs(fc(i,j,k))
               if(pcur.gt.1.)write(*,*) pcur,rk/fk
            enddo
         enddo
      enddo

      stop
      end
