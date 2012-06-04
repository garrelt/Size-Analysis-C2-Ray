typedef struct parlist {
  long seed;			/* seed for random number generator */
  float sigma;			/* width of Gaussian */
  int*dim;			/* size of grid in pixels */
  float a;			/* grid constant */
  int med;			/* number of interpolation levels */
  char*inname;			/* name of input file for data */
  char*outname;			/* name of output file for results */
  char*colname;			/* name of output file for color */
  char*grsname;			/* name of output file for greyscale */
  float lo;			/* lowest threshold */
  float hi;			/* highest threshold */
  int pixel;			/* number of pixels per cell */
  int bins;			/* number of threshold bins */
  int length;			/* number of fits to do */
  int nongauss;			/* deviation from Gaussianity */
  int time;			/* flag for giving time usage */
  int normal;			/* flag for normalizing field */
} parlist;

extern void Explain(FILE*,char*,parlist*);

