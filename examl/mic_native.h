#ifndef MIC_NATIVE_H_
#define MIC_NATIVE_H_

void printThreadCount();
void printSum(const char* label, double *arr, int n);
void printSumi(const char* label, int *arr, int n);

void mic_newviewGTRGAMMA(int tipCase,
                  double *x1, double *x2, double *x3, double *extEV, double *tipVector,
                  unsigned char *tipX1, unsigned char *tipX2,
                  int n, double *left, double *right, int *wgt, int *scalerIncrement);

double mic_evaluateGAMMA(int *wptr,
                 double *x1_start, double *x2_start,
                 double *tipVector,
                 unsigned char *tipX1, const int n, double *diagptable);

void mic_sumGAMMA(int tipCase, double *sumtable, double *x1_start, double *x2_start, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n);

void mic_coreGTRGAMMA(const int upper, double *sumtable,
    volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *EIGN, double *gammaRates, double lz, int *wrptr);

void mic_sumcoreGTRGAMMA(int tipCase, double *x1_start, double *x2_start, double *tipVector,
        unsigned char *tipX1, unsigned char *tipX2, const int n,
        volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *EIGN, double *gammaRates, double lz, int *wgt);


#endif /* MIC_NATIVE_H_ */
