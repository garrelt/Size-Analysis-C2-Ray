      program dndv

c     THIS PROGRAM CALCULATES A CATALOG OF IONIZED BUBBLES
c     USING THE EQUIVALENCE CLASS METHOD OF PRESS ET AL. CHAPTER 8
c     THE INPUT IS IONIZED FRACTION ON A MESH, THE OUTPUT IS A LIST OF
c     BUBBLES AND THEIR VOLUMES

      character *128 filein,fileout,filehalos,filegas
      parameter (nmx=406,nmx3=nmx**3,pi=3.14159)

      real, allocatable:: xarr(:,:,:),collfrac(:,:,:),rho(:,:,:)
      double precision, allocatable:: fracinclass(:),massinclass(:)
      logical, allocatable:: ionized(:),mergecell(:)
      integer, allocatable:: ijk2m(:,:,:),m2ijk(:,:),class(:),
     &         mergemap(:),nmemb(:)

      integer HIflag
      double precision tcollfrac,rhobar,xbar
      logical isolated
      real boxsize, h 
c     GET COMMANDLINE ARGUMENTS OR EXIT WITH USAGE STATEMENT

      if(iargc().le.1) then
         write(*,99)
 99      format('usage: fof_bubbles <fin> <fout> ',
     &          '[<xth> <HIfl> <Lbox> <h> <isoflag>]')
         goto 900
      endif

      call getarg(1,filein)
      call getarg(2,fileout)
      xth = r4arg(3,0.5)
      HIflag = i4arg(4,0)
      boxsize = i4arg(5,1.e2)
      h       = i4arg(6,7.e-1)
      isoflag = i4arg(7,0)
      isolated=.false.
      if(isoflag.eq.1) isolated=.true.

c     OPEN INPUT/OUTPUT FILES

      open(unit=1,file=filein,form='unformatted')
c      open(unit=1,file=filein,form='binary')
      open(unit=2,file=fileout)

c     READ INPUT

      read(1) nx,ny,nz
      nx3=nx**3

      allocate (xarr(nx,nx,nx))
      allocate (rho(nx,nx,nx))
      allocate (collfrac(nx,nx,nx))
      allocate (ionized(nx3))
      allocate (mergecell(nx3))
      allocate (ijk2m(nx,nx,nx))
      allocate (m2ijk(3,nx3))
      allocate (class(nx3))
      allocate (mergemap(nx3))
      allocate (nmemb(nx3))
      allocate (fracinclass(nx3))
      allocate (massinclass(nx3))

      read(1) xarr
c      read(3) rho

      close(1)
c      close(3)

c     (DON'T) READ SOURCE CATALOG

      collfrac=0.
c      open(unit=1,file=filehalos)
c      read(1,*) nhalos
      nhalos=0 
      tcollfrac=0.0
      gridmass=3248.**3.
      do ih=1,nhalos
c         read(1,*) i,j,k,collmass
         ii=i
         jj=j
         kk=k
         collfrac(ii,jj,kk)=collmass/gridmass
         tcollfrac=tcollfrac+collfrac(ii,jj,kk)
      enddo
c      write(*,*) 'total collapsed fraction: ',tcollfrac

      if(HIflag.eq.1) xarr=1.-xarr ! DO THIS FOR NEUTRAL REGIONS

c     INITIALIZE BOOKKEEPING AND CLASSES

      ncell=nx*ny*nz
      nion=0
      m=0
      rhobar=0.
      xbar=0.
      do i=1,nx
         do j=1,ny
            do k=1,nz
               m=m+1
               ijk2m(i,j,k)=m
               m2ijk(1,m)=i
               m2ijk(2,m)=j
               m2ijk(3,m)=k
               class(m)=m
               if(xarr(i,j,k).gt.xth) then
                  ionized(m)=.true.
                  nion=nion+1
               else
                  ionized(m)=.false.
               endif
c               rhobar=rhobar+rho(i,j,k)
c               xbar=xbar+xarr(i,j,k)*rho(i,j,k)
            enddo
         enddo
      enddo
c      rhobar=rhobar/real(ncell)
c      xbar=xbar/real(ncell)/rhobar

c      write(*,*) 'fraction cells > x_th, ionized mass fraction: ',
c     &            real(nion)/real(ncell),xbar
c      write(*,*) 'ratio of collapsed and ionized mass fractions: ',
c     &            xbar/tcollfrac

c     LOOP OVER EACH CELL'S NEIGHBORS TO DETERMINE FRIENDSHIP AND MERGE CLASSES

c      rho=rho/rhobar
      do m=1,ncell

c        GET GRID POINT OF CELL

         ii=m2ijk(1,m)
         jj=m2ijk(2,m)
         kk=m2ijk(3,m)

c        GET CLASS OF CELL

         ic1=m 
 4       if(ic1.ne.class(ic1)) then
            ic1=class(ic1)
            goto 4
         endif

         do i=ii-1,ii+1
            if(isolated.and.(i.gt.nx.or.i.lt.1)) goto 31
            ic=i
            if(i.gt.nx) ic=ic-nx
            if(i.lt.1) ic=ic+nx

            do j=jj-1,jj+1
               if(isolated.and.(j.gt.nx.or.j.lt.1)) goto 30
               jc=j
               if(j.gt.ny) jc=jc-ny
               if(j.lt.1) jc=jc+ny

               do k=kk-1,kk+1
                  if(isolated.and.(k.gt.nx.or.k.lt.1)) goto 29
                  kc=k
                  if(k.gt.nz) kc=kc-nz
                  if(k.lt.1) kc=kc+nz

c                 DO COMPARISON BETWEEN CELL AND NEIGHBOR

                  mijk=ijk2m(ic,jc,kc)
                  if(mijk.eq.m) goto 29 ! EXIT IF SAME
                  
c                 IF FRIENDS THEN MAKE NEIGHBOR CLASS SAME AS THAT OF CELL

                  idiff=iabs(ii-i)+iabs(jj-j)+iabs(kk-k)!idiff<2 ==> ONLY FACES
                  if(ionized(mijk).and.ionized(m).and.idiff.lt.2) then
                     ic2=mijk
 3                   if(ic2.ne.class(ic2)) then
                        ic2=class(ic2)
                        goto 3
                     endif
                     class(ic2)=ic1
                  endif

 29            enddo
 30         enddo
 31      enddo

      enddo

c     FINAL SWEEP UP TO MERGE CLASSES

      do m=1,ncell
 9       if(class(m).ne.class(class(m))) then
            class(m)=class(class(m))
            goto 9
         endif
      enddo

c     INITIALIZE ARRAY KEEPING TRACK OF MEMBERS OF EACH BUBBLE

      do i=1,ncell
         nmemb(i)=0
      enddo

c     ACCUMULATE NUMBERS AND COLLAPSED FRACTION IN EACH CLASS

      fracinclass=0.
      do i=1,ncell
         ii=m2ijk(1,i)
         jj=m2ijk(2,i)
         kk=m2ijk(3,i)
         fcur=rho(ii,jj,kk)/real(ncell)

c        NUMBER OF MEMBERS
         nmemb(class(i))=nmemb(class(i))+1
c        COLLAPSED FRACTION
         fracinclass(class(i))=fracinclass(class(i))+collfrac(ii,jj,kk)
c        MASS CONTAINED WITHIN
         massinclass(class(i))=massinclass(class(i))+fcur

      enddo

      nbubbles=0
      nionized=0
      boxsize=boxsize/h
      do i=1,ncell
         if(nmemb(i).gt.1) then
            nbubbles=nbubbles+1
            vol=boxsize**3*real(nmemb(i))/real(ncell)
            rad=(3.*vol/4./pi)**(1./3.)
            write(2,10) nbubbles,vol/boxsize**3,rad,fracinclass(i),
     &                  massinclass(i)
            nionized=nionized+nmemb(i)
         elseif(nmemb(i).eq.1.and.ionized(i)) then
            nbubbles=nbubbles+1
            vol=boxsize**3*real(nmemb(i))/real(ncell)
            rad=(3.*vol/4./pi)**(1./3.)
            write(2,10) nbubbles,vol/boxsize**3,rad,fracinclass(i),
     &                  massinclass(i)
            nionized=nionized+nmemb(i)
         endif
      enddo
 10   format(i7,4(1x,1pe13.6))

      close(2)

      write(*,101) filein(9:14),real(nionized)/real(ncell)
 101  format(a6,1x,f10.7)

 900  continue

      stop
      end
