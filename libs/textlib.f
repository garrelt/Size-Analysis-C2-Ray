c-----------------------------------------------------------
c     Utility routines used in Phy329 course and available 
c     in 'libp329f.a' on einstein.
c
c     Copyright Matthew W. Choptuik, 1988-1997
c
c     User routines:
c
c        i4arg
c        r4arg
c        getu
c        indlnb
c
c     Support routines
c
c        onedot
c        s2i4
c        s2r8
c-----------------------------------------------------------

c-----------------------------------------------------------
c     Converts argument 'argno' to integer*4 and returns 
c     value or 'defval' if parsing error or if argument is 
c     single dot ('.')
c-----------------------------------------------------------
      integer function i4arg(argno,defval)

         implicit       none
         
         real*8         r8_never
         parameter    ( r8_never = -1.0d-60 )

         integer        i4_never 
         parameter    ( i4_never = -2 000 000 000 )

         logical        onedot
         integer        iargc
         integer        s2i4

         integer        argno
         integer        defval

         character*32   argi

         if( argno .le. iargc() ) then
            call getarg(argno,argi)
            if( onedot(argi) ) then
               i4arg = defval
            else
               i4arg = s2i4(argi)
               if( i4arg .eq. i4_never ) then
                  i4arg = defval
               end if
            end if
         else
            i4arg = defval
         end if

         return

      end

c-----------------------------------------------------------
c     Converts argument 'argno' to real*8 and returns 
c     value or 'defval' if parsing error or if argument is 
c     single dot ('.')
c-----------------------------------------------------------
      real function r4arg(argno,defval)

         implicit       none
         
         real         r8_never
         parameter    ( r8_never = -1.0e-20 )

         integer        i4_never 
         parameter    ( i4_never = -2 000 000 000 )

         logical        onedot 
         integer        iargc
         real*8         s2r8

         integer        argno
         real         defval

         character*32   argi

         if( argno .le. iargc() ) then
            call getarg(argno,argi)
            if( onedot(argi) ) then
               r4arg = defval
            else 
               r4arg = s2r8(argi)
               if( r4arg .eq. r8_never ) then
                  r4arg = defval
               end if
            end if
         else
            r4arg = defval
         end if

         return

      end
 
c-----------------------------------------------------------

      double precision function r8arg(argno,defval)

         implicit       none

         real*8         r8_never
         parameter    ( r8_never = -1.0d-60 )

         integer        i4_never
         parameter    ( i4_never = -2 000 000 000 )

         logical        onedot
         integer        iargc
         real*8         s2r8

         integer        argno
         real*8         defval

         character*32   argi   

         if( argno .le. iargc() ) then
            call getarg(argno,argi)
            if( onedot(argi) ) then
               r8arg = defval
            else
               r8arg = s2r8(argi)
               if( r8arg .eq. r8_never ) then
                  r8arg = defval
               end if
            end if
         else
            r8arg = defval
         end if

         return

      end
                         
c-----------------------------------------------------------
c     Returns index of last non-blank character in S.
c-----------------------------------------------------------
      integer function indlnb(s)
 
         character*(*)    s
 
         do indlnb = len(s) , 1 , -1
            if( s(indlnb:indlnb) .ne. ' ' ) RETURN
         end do
         indlnb = 0
 
         return
 
      end

c-----------------------------------------------------------
c     Returns first unit number .ge. umin not attached to 
c     a file.
c-----------------------------------------------------------
      integer function getu()

         implicit      none

         integer       umin
         parameter   ( umin = 10 )

         integer       umax
         parameter   ( umax = 99 )

         integer       u
         logical       opened

         getu = -1
         do u = umin , umax
            inquire(unit=u,opened=opened)
            if( .not. opened ) then
               getu = u
               return
            end if       
         end do
         write(0,*) 'getu: Panic--no available unit number'
         stop

      end

c-----------------------------------------------------------
c     Utility routine for parameter handling.
c-----------------------------------------------------------
      logical function onedot(s)

         character*(*)     s

         integer           indlnb

         onedot = s(1:1) .eq. '.'  .and.  indlnb(s) .eq. 1

         return

      end

c-----------------------------------------------------------
c     Converts string to integer.
c-----------------------------------------------------------
      integer function s2i4(s)

         implicit      none

         real*8       r8_never
         parameter  ( r8_never = -1.0d-60 )

         integer      i4_never 
         parameter  ( i4_never = -2 000 000 000 )

         integer       indlnb

         character*(*) s     

         character*32  buffer
         integer       lens

         integer       default
         parameter   ( default = 0 )

         lens = indlnb(s)
         if( lens .gt. 99 ) then
            write(*,*) '>>> s2i4:: String too long for conversion.'
            s2i4 = default
         else 
            if( lens .le. 9 ) then
               write(buffer,100) lens
            else 
               write(buffer,101) lens
            end if
 100        format('(I',i1,')')
 101        format('(I',i2,')')
            read(s,fmt=buffer,end=900,err=900) s2i4
         end if

         return

 900        s2i4 = i4_never
         return

      end

c-----------------------------------------------------------
c     Converts string to real*8 value.
c-----------------------------------------------------------
      double precision function s2r8(s)

         implicit      none

         real*8       r8_never
         parameter  ( r8_never = -1.0d-60 )

         integer      i4_never 
         parameter  ( i4_never = -2 000 000 000 )

         integer       indlnb

         character*(*) s     

         character*32  buffer
         integer       lens

         double precision default
         parameter      ( default = 0.0d0 )

         lens = indlnb(s)
         if( lens .gt. 99 ) then
            write(*,*) '>>> s2r8:: String too long for conversion.'
            s2r8 = default
         else 
            if( lens .le. 9 ) then
               write(buffer,100) lens
            else 
               write(buffer,101) lens
            end if
 100        format('(G',i1,'.0)')
 101        format('(G',i2,'.0)')
            read(s,fmt=buffer,end=900,err=900) s2r8
         end if

         return

 900        s2r8 = r8_never 
         return

      end
