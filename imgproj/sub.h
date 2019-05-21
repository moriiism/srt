#ifndef MORIIISM_SRT_IMGPROJ_SUB_H_
#define MORIIISM_SRT_IMGPROJ_SUB_H_

#include "mir_hist1d_nerr.h"
#include "mir_hist2d_nerr.h"
#include "mir_vect.h"

void FillTriLo(double value_bin, double area_bin,
               const Vect2d* const vect_c0,
               const Vect2d* const vect_c1,
               const Vect2d* const vect_c2,
               double y_lo, double y_up,
               HistDataNerr1d* const hd1d_proj);

void FillTriUp(double value_bin, double area_bin,
               const Vect2d* const vect_c1,
               const Vect2d* const vect_c2,
               const Vect2d* const vect_c3,
               double y_lo, double y_up,
               HistDataNerr1d* const hd1d_proj);

void FillSq(double value_bin, double area_bin,
            const Vect2d* const vect_c0,
            const Vect2d* const vect_c1,
            const Vect2d* const vect_c2,
            const Vect2d* const vect_c3,
            double y_lo, double y_up,
            HistDataNerr1d* const hd1d_proj);

double GetCrossYPos(double xpos1, double ypos1,
                    double xpos2, double ypos2,
                    double xpos);

#endif // MORIIISM_SRT_IMGPROJ_SUB_H_
