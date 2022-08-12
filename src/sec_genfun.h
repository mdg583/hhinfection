#ifndef SEC_GENFUN_H
#define SEC_GENFUN_H

#include <math.h>

struct Gdat {
    double q;
    double *data;
    int max_hh;
};
Gdat gdat(int max_hh, double q);
Gdat lgdat(int max_hh, double q);
long gsize(Gdat dat);
double *g(Gdat gdat, int k, int m);
void gfree(Gdat gdat);

#endif /* SEC_GENFUN_H */
