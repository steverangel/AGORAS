#ifndef _H_CLARANS
#define _H_CLARANS
float searchNeighborhood(int , int, int,float **, int *, int *,float *, int*,float(*)(int,float *, float *), int, float**);
/* called from main() */
int *seq_clarans(float **, int, int, int, int *, float(*distcalc)(int,float *, float *),int, int, float*, unsigned int);
#endif
