#ifndef MORIIISM_SRT_MKBG_MODEL_H_
#define MORIIISM_SRT_MKBG_MODEL_H_

#include "mir_math.h"

double GetConstFunc(double posx, double posy, double constant);
double GetLinFunc(double posx, double posy, double w0, double w1, double constant);
double GetQuadFunc(double posx, double posy,
                   double posx_ctr, double posy_ctr,
                   double coeff0, double coeff1, double coeff2,
                   double constant);

#endif // MORIIISM_SRT_MKBG_MODEL_H_
