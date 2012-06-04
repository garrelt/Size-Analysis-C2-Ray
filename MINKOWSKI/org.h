extern void empty(char**);
extern void fileopenr(FILE**,char*);
extern void fileopenw(FILE**,char*);
extern void fileclose(FILE*);

extern double chisquare(int,long*);

extern int*intvector(int,int);
extern float*vector(int,int);
extern double*doublevector(int,int);
extern int**intmatrix(int,int,int,int);
extern float**matrix(int,int,int,int);
extern double**doublematrix(int,int,int,int);
extern float***cube(int,int,int,int,int,int);
extern void free_intvector(int*,int,int);
extern void free_vector(float*,int,int);
extern void free_doublevector(double*,int,int);
extern void free_intmatrix(int**,int,int,int,int);
extern void free_doublematrix(double**,int,int,int,int);
extern void free_matrix(float**,int,int,int,int);
extern double gasdev(long*);
extern double ran4(long*);
extern float qromb(float(*)(float),float,float);
