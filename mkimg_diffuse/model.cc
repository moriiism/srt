#include "model.h"

double GetConstFunc(double posx, double posy, double constant)
{
    double val = constant;
    return(val);
}

double GetLinFunc(double posx, double posy,
                  double w0, double w1, double constant)
{
    double val = 0.0;
    val = w0 * posx + w1 * posy + constant;
    return(val);
}

double GetQuadFunc(double posx, double posy,
                   double posx_ctr, double posy_ctr,
                   double coeff0, double coeff1, double coeff2,
                   double constant)
{
    double val = coeff0 * (pow(posx - posx_ctr, 2) + pow(posy - posy_ctr, 2))
        + coeff1 * posx + coeff2 * posy + constant;
    return(val);
}
