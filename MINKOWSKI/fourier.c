#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "main.h"
#include "minkowski.h"
#include "org.h"
#define pi 3.141592653589793238462643

void convolution(float*,parlist*);
void randomfield(float*,parlist*);
void gaussiankernel(float*,parlist*);
void gaussian(float*,int,float);
void normalize(float*,parlist*);

void fourn(float*,int*,int,int);
int idx(int,int,int,int*);

/* --------------------------------------------- */
/*                                               */
/* routines for generating fields and kernels    */
/* --------------------------------------------- */

/* Random field with scale free Power spectrum */
void randomfield(float*u,parlist*par)
{
  int i,j,k,l, ii,jj,kk;
  int lx=par->dim[1], ly=par->dim[2], lz=par->dim[3];
  float real,imag,sigma,ran, x,y,z;

  /* make a random field */
  for(i=0; i<2*lx*ly*lz; i++) u[i]=0;
  for(i=0; i<=lx/2; i++) for(j=0; j<ly; j++) for(k=0; k<lz; k++)
    if(i<lx/2 && i>0 || j<ly/2) {
      ii = i?lx-i:0; x=i/(float)lx;
      jj = j?ly-j:0; y=(j<jj?j:jj)/(float)ly;
      kk = k?lz-k:0; z=(k<kk?k:kk)/(float)lz;
      if(x||y||z) {
	sigma=sqrt(x*x+y*y+z*z);
	if(par->nongauss) {
	  real=sigma*chisquare(par->nongauss,&par->seed);
	  imag=sigma*chisquare(par->nongauss,&par->seed); 
	} else {
	  real=sigma*gasdev(&par->seed);
	  imag=sigma*gasdev(&par->seed);
	}
	if(i==lx/2||j==ly/2||k==lz/2) imag=0;
      }
      else real=imag=0;
      u[idx( i, j, k,par->dim)]=+real, u[idx( i, j, k,par->dim)+1]=+imag;
      u[idx(ii,jj,kk,par->dim)]=+real, u[idx(ii,jj,kk,par->dim)+1]=-imag;
    }
}

/* convolution with Gaussian kernel */
void convolution(float*u,parlist*par)
{
  int i,j,k,*d=par->dim;
  float*x,*y,*z;

  /* allocate memory and calculate one-dimensional Gaussians */
  x=vector(0,d[1]); gaussian(x,d[1],par->sigma*2*pi/d[1]);
  y=vector(0,d[2]); gaussian(y,d[2],par->sigma*2*pi/d[2]);
  z=vector(0,d[3]); gaussian(z,d[3],par->sigma*2*pi/d[3]);

  /* convolution in Fourier space */
  for(i=0;i<d[1];i++) for(j=0;j<d[2];j++) for(k=0;k<d[3];k++)
    u[idx(i,j,k,d)  ]*=x[i]*y[j]*z[k],
    u[idx(i,j,k,d)+1]*=x[i]*y[j]*z[k];
} 

/* Gaussian kernel in real space */
/* sigma -- well, guess what     */
void gaussiankernel(float*g, parlist*par)
{
  int i,j,k,*d=par->dim;
  float sig=par->sigma, norm, *x,*y,*z;

  /* allocate memory and calculate one-dimensional Gaussians */
  x=vector(0,d[1]); gaussian(x,d[1],1/sig);
  y=vector(0,d[2]); gaussian(y,d[2],1/sig);
  z=vector(0,d[3]); gaussian(z,d[3],1/sig);

  /* calculate multi-dimensional Gaussian */
  norm=pow(2*pi*sig*sig,1.5);
  for(i=0;i<d[1];i++) for(j=0;j<d[2];j++) for(k=0;k<d[3];k++)
    g[idx(i,j,k,d)]=x[i]*y[j]*z[k]/norm, g[idx(i,j,k,d)+1]=0;
  
  /* free memory */
  free_vector(x,0,d[1]);
  free_vector(y,0,d[2]);
  free_vector(z,0,d[3]);
}
/* Gaussian in one dimension */
void gaussian(float*x,int n,float sig)
{
  int i,ii;

  for(i=0; i<n; i++) {
    if(i>n/2) ii=i-n; else ii=i;
    x[i]=exp(-.5*ii*ii*sig*sig);
  }

  return;
}

/* normalize field */
void normalize(float*u,parlist*par)
{
  int i,j,k,*d=par->dim,n=d[1]*d[2]*d[3];
  double dmu=0,dsig=0;
  float mu=0,sig=1,value;
  float min=1.e10,max=-min;

  /* normalize to <f>=0, <f^2>=1 */
  if(par->normal) {
    mu=sig=0;
    for(i=0;i<d[1];i++) for(j=0;j<d[2];j++) for(k=0;k<d[3];k++)
      /*  value=u[idx(i,j,k,d)], mu+=value,sig+=value*value;  /*was commented out*/
      value=u[idx(i,j,k,d)], dmu+=value,dsig+=value*value; 
    mu=dmu;sig=dsig;
    sumtomean(&mu,&sig,n);
    for(i=0;i<n;i++){
      if(u[2*i]<min) min=u[2*i];
      if(u[2*i]>max) max=u[2*i];
    }
/*   for(i=0;i<n;i++) u[2*i]=(u[2*i]-mu)/sig; /* NOT NORMALIZING!!! */ 
 /*   printf("mu,sig,min,max = %g %g %g %g\n",mu,sig,min,max);   */
  }

}

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void fourn(float*data,int*nn,int ndim,int isign)
{
  int i;
  int i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
  int ibit,idim,k1,k2,n,nprev,nrem,ntot;
  float tempi,tempr;
  float theta,wi,wpi,wpr,wr,wtemp;

  ntot=1;
  for (idim=1;idim<=ndim;idim++) ntot *= nn[idim];
  nprev=1;
  for (idim=ndim;idim>=1;idim--) {
    n=nn[idim];
    nrem=ntot/(n*nprev);
    ip1=nprev << 1;
    ip2=ip1*n;
    ip3=ip2*nrem;
    i2rev=1;
    for (i2=1;i2<=ip2;i2+=ip1) {
      if (i2 < i2rev) {
	for (i1=i2;i1<=i2+ip1-2;i1+=2) {
	  for (i3=i1;i3<=ip3;i3+=ip2) {
	    i3rev=i2rev+i3-i2;
	    SWAP(data[i3],data[i3rev]);
	    SWAP(data[i3+1],data[i3rev+1]);
	  }
	}
      }
      ibit=ip2 >> 1;
      while (ibit >= ip1 && i2rev > ibit) {
	i2rev -= ibit;
	ibit >>= 1;
      }
      i2rev += ibit;
    }
    ifp1=ip1;
    while (ifp1 < ip2) {
      ifp2=ifp1 << 1;
      theta=isign*6.28318530717959/(ifp2/ip1);
      wtemp=sin(0.5*theta);
      wpr=-2.0*wtemp*wtemp;
      wpi=sin(theta);
      wr=1.0;
      wi=0.0;
      for (i3=1;i3<=ifp1;i3+=ip1) {
	for (i1=i3;i1<=i3+ip1-2;i1+=2) {
	  for (i2=i1;i2<=ip3;i2+=ifp2) {   /*added -ifp1*/
	    k1=i2;
	    k2=k1+ifp1;
	    tempr=wr*data[k2]-wi*data[k2+1];
	    tempi=wr*data[k2+1]+wi*data[k2];
	    data[k2]=data[k1]-tempr;
	    data[k2+1]=data[k1+1]-tempi;
	    data[k1] += tempr;
	    data[k1+1] += tempi;
	  }
	}
	wr=(wtemp=wr)*wpr-wi*wpi+wr;
	wi=wi*wpr+wtemp*wpi+wi;
      }
      ifp1=ifp2;
    }
    nprev *= n;
  }

  /* normalize the inverse transform */
  if(isign==-1) for(i=1;i<=2*ntot;i++) data[i]/=ntot;
}

#undef SWAP

int idx(int i,int j,int k,int*dim)
{
  while(i>=dim[1]) i-=dim[1]; while(i<0) i+=dim[1];
  while(j>=dim[2]) j-=dim[2]; while(j<0) j+=dim[2];
  while(k>=dim[3]) k-=dim[3]; while(k<0) k+=dim[3];
  return 2*((i*dim[2]+j)*dim[3]+k);
}
