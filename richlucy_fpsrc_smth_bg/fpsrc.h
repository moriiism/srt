#ifndef MORIIISM_SRT_RICHLUCY_FPSRC_SMTH_BG_FPSRC_H_
#define MORIIISM_SRT_RICHLUCY_FPSRC_SMTH_BG_FPSRC_H_

#include "mi_iolib.h"
#include "mib_blas.h"
#include "mir_math.h"

void GenFixedPointSrcDetImg(string fixed_src_list,
                            const double* const resp_norm_mat_arr,
                            int nskyx, int nskyy, int ndet,
                            int* const nsrc_ptr,
                            double*** const det_fpsrc_arr_ptr);

#endif // MORIIISM_SRT_RICHLUCY_FPSRC_SMTH_BG_FPSRC_H_
