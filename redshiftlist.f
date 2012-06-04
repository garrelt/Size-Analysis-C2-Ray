      program redshiftlist

c     THIS PROGRAM LOADS THE IONIZATION FRACTION FILES 
c     AND CALCULATES THE VOLUME AND MASS AVERAGED IONIZATION
c     FRACTION

      character *256 filein1,filein2,fileout

      real, allocatable:: x(:,:,:)
      real, allocatable:: ndens(:,:,:)
      integer:: n, i, j, k 
      real :: out(2)

      if(iargc().le.2) then
         write(*,99)
 99      format('usage: redshiftlist <fin1> <dens in> <fout>')
         goto 900
      endif
           
      call getarg(1,filein1)
      call getarg(2,filein2)
      call getarg(3,fileout)
      n=256
      open(unit=1,file=filein1,form='unformatted')
      read(1) n,n,n

      allocate (x(n,n,n))
      allocate (ndens(n,n,n))


      read(1) (((x(i,j,k),k=1,n),j=1,n),i=1,n)
      close(1)    

      open(unit=1,file=filein2,form='binary')
      read(1) n,n,n
      read(1) (((ndens(i,j,k),k=1,n),j=1,n),i=1,n)
      close(1)

      xtot=0.
      xtot2=0.
      dentot=0.
      do i=1,n
      do j=1,n
      do k=1,n
         xtot=xtot+x(i,j,k)
         xtot2=xtot2+x(i,j,k)*ndens(i,j,k)
         dentot=dentot+ndens(i,j,k)
      enddo
      enddo
      enddo
         xtot=xtot/(n*n*n)
         xtot2=xtot2/dentot   
         
         out(1)=xtot
         out(2)=xtot2

      open(unit=52,file=fileout,status="unknown")
      write(52,11) (out(k),k=1,2)
 11   format(2(1pe13.6,1x))
      close(52)

 900  continue

      stop
      end

