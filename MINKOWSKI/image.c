#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "main.h"
#include "org.h"
#define pi 3.141592653589793238462643
#define NPOINTS 20000
#define true 1
#define false 0

void picture(int, FILE*, float*, parlist*);
void bitmapxpm(int, FILE*, char**, int,int);
char**cmatrix(int,int);

void picture(int cl, FILE*fp, float*z, parlist*par)
{
  int i,j,k, ii,jj, ip,jp, c,r, free;
  int *d=par->dim,n=par->dim[1]*par->dim[2]*par->dim[3];
  int nx,ny, lx,ly, step,pixl;
  char**map;
  float zmin,zmax, alpha;

  /* size of cell in pixels */
  if(par->pixel>0) pixl=par->pixel, step=1; else pixl=1, step=-par->pixel;
  
  /* calculate range of points */
  zmin=66666; zmax=-66666;
  for(i=0;i<n;i++) {
    if(z[2*i]!=-666&&z[2*i]>zmax) zmax=z[2*i];
    if(z[2*i]!=-666&&z[2*i]<zmin) zmin=z[2*i];
  }
/*  fprintf(stderr,"%d data points range from %g to %g\n",n,zmin,zmax);*/
  zmin=par->lo; zmax=par->hi;

  /* initialize bitmap */
  for(nx=sqrt(d[3]/step); d[3]/step%nx; nx--); ny=d[3]/step/nx;
  lx=nx*(d[1]*pixl/step+10)+10, ly=ny*(d[2]*pixl/step+10)+10; 
  map=cmatrix(lx,ly); for(i=0;i<lx;i++) for(j=0;j<ly;j++) map[i][j]=48;

  /* put points in matrix */
  for(i=0;i<d[1];i+=step) for(j=0;j<d[2];j+=step) for(k=0;k<d[3];k+=step) {
    c=63*(z[idx(i,j,k,d)]-zmin)/(zmax-zmin); if(c<0) c=0; if(c>63) c=63;
    ii=i/step+(k/step%nx)*(d[1]*pixl/step+10)+10; 
    jj=j/step+(k/step/nx)*(d[2]*pixl/step+10)+10;
    for(ip=0;ip<pixl;ip++) for(jp=0;jp<pixl;jp++) map[ii+ip][jj+jp]=c%64+49;
  } 

  /* output bitmap */
  bitmapxpm(cl,fp,map,lx,ly);
}

/* output matrix in XPM bitmap format */
void bitmapxpm(int cl,FILE*fp,char**bitmap,int lx,int ly)
{
  int i,j, x,y, r;

  /* header and colour table */
  fprintf(fp,"/* XPM */\n"
	  "static char * bitmap [] = {\n"
	  "/* size height ncolors cpp [x_hot y_hot] */\n"
	  "\"%d %d 65 1 0 0\",\n"
	  "/* colours */\n"
	  "\"%c\tc #ffffff\",\n", 
	  /* "\"%c\tc #000000\",\n",  */
	  lx,ly,48);
  for(i=49; i<=112; i++) {
    fprintf(fp,"\"%c\tc #",i);
    if(cl) {
      /* blue(cold) to red(hot) */
      if(i<=52)fprintf(fp,"0"); fprintf(fp,"%x",4*(i-49)+3);        /* R */
      fprintf(fp,"00");                                             /* G */
      if(i>=109)fprintf(fp,"0"); fprintf(fp,"%x\",\n",4*(112-i)+3); /* B */
    } else {
      /* greyscale white(cold) to black(hot) */
      for(j=0; j<3; j++) {
	if(i>=109) fprintf(fp,"0"); fprintf(fp,"%x",4*(112-i)+3);
      }
      fprintf(fp,"\",\n");
    }
  }

  /* pixels */
  fprintf(fp,"/* pixels */\n");
  for(i=0;i<ly;i++) {
    fprintf(fp,"\"");
    for(j=0;j<lx;j++) fprintf(fp,"%c",bitmap[j][i]);
    if(i==ly-1) fprintf(fp,"\"};\n"); else fprintf(fp,"\",\n");
  }
}

/* allocate memory for matrix of chars */
char**cmatrix(int width, int height)
{
  int i;
  char**u;

  u=(char**)malloc(width*sizeof(char*));
  for(i=0; i<width; i++) u[i]=(char*)calloc(height,sizeof(char));

  return u;
}
