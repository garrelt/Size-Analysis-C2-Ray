cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine FFT3RKs(fc,n1,n2,n3)
c  From real space to k-space on a single processor.
c  fc = input array of reals packed into half-size complex with one extra
c       x-plane (the Nyquist plane).  fc is returned as the desired transform.
c  (n1,n2,n3) = dimension of full real array.
c
	implicit none
	integer, intent(in):: n1,n2,n3
	complex, intent(inout)::fc(n1/2+1,n2,n3)

        integer, parameter:: isign=-1
	complex, dimension(n2,n3):: fe,fo
	complex:: z
	complex:: im
	integer:: k1,k1c,k2,k3,n12,n14,n21,j,i1
	real:: theta
	real, parameter:: twopi=6.283185307179586d0

	im=cmplx(0.0,1.0)
c
	n12=n1/2
	n14=n1/4
	n21=n12+1

c  FFT 1st dimension of a 3-D array.  Note that first dimension has
c  length n12+1 but transform length is only n12.
	call FFT13(fc,n12,n2,n3,isign)

c  Unpack transform of real data in the 1st dimension.
c  First, zero and Nyquist frequencies.
c  Pack fcny and fc(1,:,:) into fc(1,:,:) so we can do remaining ffts
c  in one shot.
	fc(1,:,:)=cmplx(real(fc(1,:,:))+aimag(fc(1,:,:)),
     &                 real(fc(1,:,:))-aimag(fc(1,:,:)))

c  Unpack other planes.
	do k1=2,n14+1
	  k1c=n12+2-k1
	  theta=isign*twopi*(k1-1)/n1
	  z=cmplx(cos(theta),sin(theta))
	  fe(:,:)=0.5*(fc(k1,:,:)+conjg(fc(k1c,:,:)))
	  fo(:,:)=0.5*(fc(k1,:,:)-conjg(fc(k1c,:,:)))
	  fo(:,:)=-im*z*fo(:,:)
	  fc(k1,:,:)=fe(:,:)+fo(:,:)
	  fc(k1c,:,:)=conjg(fe(:,:)-fo(:,:))
	end do

c  Do 2-D FFTs on 2nd and 3rd dimensions.
c  First transform 2nd dimension.
	call FFT23(fc,n21,n2,n3,isign)
c  Next transform 3rd dimension.
	call FFT33(fc,n21,n2,n3,isign)

c  Now need to unscramble fc(1,:,:) into fc(1,:,:) and fc(n12+1,:,:).
	fe = fc(1,:,:)
	fo = 0.0

	fo(1,1)=fc(1,1,1)
	do k2=2,n2
	  fo(k2,1)=fc(1,n2+2-k2,1)
	end do

	do k3=2,n3
	  fo(1,k3) = fc(1,1,n3+2-k3)
	end do

	do k3=2,n3
	  do k2=2,n2
	    fo(k2,k3) = fc(1,n2+2-k2,n3+2-k3)
	  end do
	end do

	fo = conjg(fo)
	fc(1,:,:) = 0.5*(fe + fo)
	fc(n21,:,:) =  -0.5*im*(fe - fo)

        return
	end subroutine FFT3RKs
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine FFT3KRs(fc,n1,n2,n3)
c  From k-space to real space on a single processor.
c  fc = input array of complex transform values including the Nyquist plane.
c       fc is returned as the desired transform, packed into half-size complex.
c  (n1,n2,n3) = dimension of full real array.
c
	implicit none
	integer, intent(in):: n1,n2,n3
	complex, intent(inout)::fc(n1/2+1,n2,n3)

        integer, parameter:: isign=+1
	complex, dimension(n2,n3):: fe,fo
	complex:: z
	complex:: im
	integer:: k1,k1c,k2,k3,n12,n14,n21,j,i1
	real:: theta
	real, parameter:: twopi=6.283185307179586d0

	im=cmplx(0.0,1.0)
c
	n12=n1/2
	n14=n1/4
	n21=n12+1

c  Pack Nyquist plane into fc(1,:,:) so that we can FFT it along with
c  all the other planes.
	fc(1,:,:)=fc(1,:,:)+im*fc(n21,:,:)

c  Do 2-D FFTs on 2nd and 3rd dimensions.

c  First FFT 3rd dimension of fc.
	call FFT33(fc,n21,n2,n3,isign)
c  Now transform 2nd dimension.
	call FFT23(fc,n21,n2,n3,isign)

c  Pack inverse transform of real data in the 1st dimension.
c  Multiply by two so that forward and inverse transform requires
c  dividing by n1*n2*n3 as for a full 3-D complex FFT.
c  First, zero frequency.
c  fc(1,:,:) and fc(n21,:,:) are both real, but they have been packed into
c     the real and imag parts of fc(1,:,:).  So first we separate them out.
c     (We could combine this with the step below, but this makes it more clear.)
	fc(n21,:,:)=-im*0.5*(fc(1,:,:)-conjg(fc(1,:,:)))
	fc(1,:,:)=0.5*(fc(1,:,:)+conjg(fc(1,:,:)))
c  Now pack other planes.
	fe(:,:)=fc(1,:,:)+conjg(fc(n21,:,:))
	fo(:,:)=fc(1,:,:)-conjg(fc(n21,:,:))
	fo(:,:)= im*fo(:,:)
	fc(1,:,:)=fe(:,:)+fo(:,:)
	fc(n21,:,:)=cmplx(0.0,0.0)
	do  k1=2,n14+1
	   k1c=n12+2-k1
	   theta=isign*twopi*(k1-1)/n1
	   z=cmplx(cos(theta),sin(theta))
	   fe(:,:)=fc(k1,:,:)+conjg(fc(k1c,:,:))
	   fo(:,:)=fc(k1,:,:)-conjg(fc(k1c,:,:))
	   fo(:,:)=im*z*fo(:,:)
	   fc(k1,:,:)=fe(:,:)+fo(:,:)
	   fc(k1c,:,:)=conjg(fe(:,:)-fo(:,:))
	enddo

c  FFT 1st dimension of a 3-D array.  Note that first dimension has
c  length n12+1 but transform length is only n12.
	call FFT13(fc,n12,n2,n3,isign)

	return
	end subroutine FFT3KRs
