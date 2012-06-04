#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "main.h"
#include "fourier.h"
#include "minkowski.h"
#include "image.h"
#include "org.h"
#define pi 3.141592653589793238462643
#define READINT atoi(ar[n][2] ? &ar[n][2] : ar[++n])
#define READFLT atof(ar[n][2] ? &ar[n][2] : ar[++n])
#define READSTR     (ar[n][2] ? &ar[n][2] : ar[++n])

void main(int nar, char* ar[])
{
  int i,j,k,l,m, n=0, num, nnn[4];
  float*u,*f;
  parlist par;
  FILE *in,*out,*col,*grs;
  char name[100];
  float cpu[2]; 
  clock_t t;

  /* set default values */
  par.sigma=2; par.seed=200294;
  par.dim=intvector(1,3); par.dim[1]=64; par.dim[2]=32; par.dim[3]=16; 
  empty(&par.inname); empty(&par.outname);
  empty(&par.colname); empty(&par.grsname);
  par.length=1; par.nongauss=0; par.lo=-4; par.hi=4; par.bins=100; 
  par.med=0; par.pixel=1; par.normal=1; par.time=0;
  
  /* work through arguments */
  while(++n<nar){
    if(ar[n][0]=='-'&&ar[n][1]){
      switch(ar[n][1]) {
	case'0' : par.grsname =READSTR; break;
	case'1' : par.colname =READSTR; break;
	case'L' : par.length  =READINT; break;
	case'N' : par.normal  =0;       break;
	case'b' : par.bins    =READINT; break;
	case'g' : par.nongauss=READINT; break;
	case'h' : par.hi      =READFLT; break;
	case'i' : par.inname  =READSTR; break;
	case'l' : par.lo      =READFLT; break;
	case'm' : par.med     =READINT; break;
	case'o' : par.outname =READSTR; break;
	case'p' : par.pixel   =READINT; break;
	case'r' : par.seed    =READINT; break;
	case's' : par.sigma   =READFLT; break;
	case't' : par.time    =1;       break;
	case'x' : par.dim[1]  =READINT; break;
	case'y' : par.dim[2]  =READINT; break;
	case'z' : par.dim[3]  =READINT; break;
      default : Explain(stderr,ar[0],&par);
      }
    }
    else {
      Explain(stderr,ar[0],&par);
    }
  }
  par.a=1./(float)par.dim[1];	/* grid constant is used everywhere */

  /* allocate some memory */
  n=par.dim[1]*par.dim[2]*par.dim[3];
  u=vector(0,2*n);

  /* open input file and read data */
  fileopenr(&in,par.inname);
 /*   printf("%s"," here 1a \n");

    printf("size = %d\n",sizeof(in));   */
  if(in) {
    f=vector(0,n);

    if(n!=fread((void*)f,sizeof(float),n,in)) {
      fprintf(stderr,"error when reading input!\n"); 
      exit(99);
    }  
/*      printf("f[0] %d\n",&f[0]);
      printf("f[0]wert %g\n", f[0]);
      printf("f[1]wert %g\n", f[1]);
      printf("f[1]wert %g\n", f[2]);
      printf("f[1]wert %g\n", f[3]);  */
      printf("---\n");    
 
    for(i=0;i<n;i++) u[2*i]=f[i],u[2*i+1]=0; 
 /*     printf("u[0] %d\n",&u[0]);
      printf("u[0]wert %g\n", u[0]);
      printf("f[0] %d\n",&f[0]);
      printf("f[0]wert %g\n", f[0]);
      printf("u[1] %d\n",&u[1]);
      printf("u[1]wert %g\n", u[1]);
      printf("f[1] %d\n",&f[1]);
      printf("f[1]wert %g\n", f[1]);
      printf("%s"," here 1c \n");    */

    free_vector(f,0,n);
 /*     printf("%s"," here 1d \n"); 
      printf("par.sigma = %g\n",par.sigma);   */

  if(par.sigma>0) fourn(u-1,par.dim,3,1); 
  }
   /*    printf("%s"," here 1e \n"); */
  fileclose(in);
  /*   printf("%s"," here 1 \n");   */
  /* open output files */
  fileopenw(&out,par.outname);
 /*   printf("%s"," here 2 \n");   */
  cpu[0]=cpu[1]=0;
  for(num=0; num<par.length; num++) {
      
    /* random field in Fourier space */
    if(par.time) t=clock();
    if(!in) randomfield(u,&par);
    if(par.time) cpu[0]+=(clock()-t)/(float)CLOCKS_PER_SEC;
    /* convolution and normalization */
    if(par.time) t=clock();
    if(par.sigma>0) convolution(u,&par);
    if(par.sigma>0) fourn(u-1,par.dim,3,-1);
    normalize(u,&par);
    if(par.time) cpu[0]+=(clock()-t)/(float)CLOCKS_PER_SEC;
  /*   printf("%s"," here 2c \n");*/
    /* perform statistics */
    if(par.time) t=clock();
    minkowski(out,u,&par);
    if(par.time) cpu[1]+=(clock()-t)/(float)CLOCKS_PER_SEC;
  }
/*      printf("%s"," here 3 \n"); */
  if(par.time) 
    fprintf(stderr,"CPU: %13s%13s\n"
	    "      %8.2f sec %8.2f sec\n",
	    "fields","minkowski",cpu[0],cpu[1]);

  /* output xpm bitmap data */
  fileopenw(&col,par.colname); if(col) picture(1,col,u,&par); fileclose(col);
  fileopenw(&grs,par.grsname); if(grs) picture(0,grs,u,&par); fileclose(grs);

  /* finish */
  fileclose(out);
  free_vector(u,0,2*n);
  exit(0);
}

/* explain options */
void Explain(FILE*fp, char*name, parlist*act)
{
  fprintf(fp,"# %s options:\n"
	  "# \t-s -- width of Gaussian                      (actual %g)\n"
	  "# \t-r -- seed for random numbers                (actual %d)\n"
	  "# \t-x -- width  of grid in cells                (actual %d)\n"
	  "# \t-y -- height of grid in cells                (actual %d)\n"
	  "# \t-z -- depth  of grid in cells                (actual %d)\n"
	  "# \t-m -- number of interpolation levels         (actual %d)\n"
	  "# \t-i -- input file for data (binary)           (actual %s)\n"
	  "# \t-o -- output file for results                (actual %s)\n"
	  "# \t-0 -- output file for greyscale image        (actual %s)\n"
	  "# \t-1 -- output file for colour image           (actual %s)\n"
	  "# \t-L -- number of fields to simulate           (actual %d)\n"
	  "# \t-g -- deviation from Gaussianity             (actual %d)\n"
	  "# \t-l -- lowest density threshold               (actual %g)\n"
	  "# \t-h -- highest density threshold              (actual %g)\n"
	  "# \t-b -- number of threshold bins               (actual %d)\n"
	  "# \t-p -- pixels per cell, negative for thinning (actual %d)\n"
	  "# \t-N -- normalize field to <u>=0, <u^2>=1      (actual %s)\n"
	  "# \t-t -- output CPU time usage                  (actual %s)\n",
	  name, act->sigma, act->seed,
	  act->dim[1], act->dim[2], act->dim[3], act->med,
	  act->inname, act->outname, act->grsname, act->colname,
	  act->length, act->nongauss, act->lo,act->hi,act->bins,act->pixel,
	  act->normal?"yes":"no", act->time?"yes":"no");
  if(fp==stderr) exit(0);
}
