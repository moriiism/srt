#include "fpsrc.h"

void GenFixedPointSrcDetImg(FILE* const fp_log,
                            string fixed_src_list,
                            const double* const resp_norm_mat_arr,
                            int nskyx, int nskyy, int ndet,
                            int* const nsrc_ptr,
                            int** const xpos_src_arr_ptr,
                            int** const ypos_src_arr_ptr,
                            double*** const det_fpsrc_arr_ptr)
{
    int nsky = nskyx * nskyy;
    
    // fixed point source list
    string* lines_arr = NULL;
    long nline = 0;
    MiIolib::GenReadFileSkipComment(fixed_src_list,
                                    &lines_arr, &nline);
    int nsrc = nline;
    int* xpos_arr = new int [nsrc];
    int* ypos_arr = new int [nsrc];
    for(long iline = 0; iline < nline; iline ++){
        string* split_arr = NULL;
        int nsplit = 0;
        MiStr::GenSplit(lines_arr[iline], &nsplit, &split_arr);
        xpos_arr[iline] = atoi(split_arr[0].c_str());
        ypos_arr[iline] = atoi(split_arr[1].c_str());
        MiIolib::Printf2(fp_log, 
                         "%d  %d\n", xpos_arr[iline], ypos_arr[iline]);
        delete [] split_arr;
    }
    delete [] lines_arr;

    double** det_fpsrc_arr = new double*[nsrc];
    for(int isrc = 0; isrc < nsrc; isrc++){
        det_fpsrc_arr[isrc] = new double[ndet];
    }

    for(int isrc = 0; isrc < nsrc; isrc++){
        // point srcfile
        double* sky_arr = new double [nsky];
        for(int isky = 0; isky < nsky; isky ++){
            sky_arr[isky] = 0.0;
        }
        int isky_pos = ypos_arr[isrc] * nskyx + xpos_arr[isrc];
        sky_arr[isky_pos] = 1.0;

        // det_arr = R_mat %*% src_norm_arr
        char transa[1];
        strcpy(transa, "N");
        dgemv_(transa, ndet, nsky, 1.0,
               const_cast<double*>(resp_norm_mat_arr), ndet,
               const_cast<double*>(sky_arr), 1,
               0.0, det_fpsrc_arr[isrc], 1);
        delete [] sky_arr;
    }

    *nsrc_ptr = nsrc;
    *xpos_src_arr_ptr = xpos_arr;
    *ypos_src_arr_ptr = ypos_arr;
    *det_fpsrc_arr_ptr = det_fpsrc_arr;
}
