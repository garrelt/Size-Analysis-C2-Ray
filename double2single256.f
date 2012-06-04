      program double2single

c     THIS PROGRAM CONVERTS DOUBLE TO SINGLE PRECISION

      character *256 filein,fileout

      double precision, allocatable:: x(:,:,:)
      integer:: n,i,j,k

      if(iargc().le.1) then
         write(*,99)
 99      format('usage: double2single <fin> <fout>')
         goto 900
      endif
           
      call getarg(1,filein)
      call getarg(2,fileout)
      n=256
      
      allocate (x(n,n,n))

c     GET INPUT REAL DATA

      open(unit=1,file=filein,form='unformatted')
      read(1) n,n,n
      read(1) (((x(i,j,k),k=1,n),j=1,n),i=1,n)
c      close(1)
c      do i=1,n
c      do j=1,n
c      do k=1,n
c         if (x(i,j,k).gt.0.9) x(i,j,k)=1.0
c         if (x(i,j,k).le.0.01) x(i,j,k)=0.0
c      enddo
c      enddo
c      enddo

      open(unit=52,file=fileout,form="unformatted",status="unknown")
      write(52) n,n,n
      write(52) (((real(x(i,j,k)),k=1,n),j=1,n),i=1,n)
      close(52)

 900  continue

      stop
      end

