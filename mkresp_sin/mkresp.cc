#include "mir_math.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "arg_mkresp.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    ArgValMkresp* argval = new ArgValMkresp;
    argval->Init(argc, argv);
    argval->Print(stdout);

    char logfile[kLineSize];
    if( MiIolib::TestFileExist(argval->GetOutdir()) ){
        char cmd[kLineSize];
        sprintf(cmd, "mkdir -p %s", argval->GetOutdir().c_str());
        system(cmd);
    }
    sprintf(logfile, "%s/%s_%s.log",
            argval->GetOutdir().c_str(),
            argval->GetOutfileHead().c_str(),
            argval->GetProgname().c_str());
    FILE* fp_log = fopen(logfile, "w");
    MiIolib::Printf2(fp_log, "-----------------------------\n");
    argval->Print(fp_log);

    
    int nskyx = 60;
    int nskyy = 60;
    int ndetx = 60;
    int ndety = 60;
    int nsky  = nskyx * nskyy;
    int ndet  = ndetx * ndety;
   for(int iskyy = 0; iskyy < nskyy; iskyy ++){
        for(int iskyx = 0; iskyx < nskyx; iskyx ++){
            char outfile[kLineSize];
            sprintf(outfile, "%s/gimage_%3.3d_%3.3d.img", argval->GetOutdir().c_str(), iskyx, iskyy);
            int isky = nskyx * iskyy + iskyx;

            double* resp_arr = new double [ndet];
            int iresp = 0;
            for(int idet = 0; idet < ndet / 2; idet ++){
                resp_arr[iresp] = cos(2 * M_PI * idet * isky / nsky);
                iresp ++;
            }
            for(int idet = 1; idet <= ndet / 2; idet ++){
                resp_arr[iresp] = sin(2 * M_PI * idet * isky / nsky);
                iresp ++;
            }
            for(int idet = 0; idet < ndet; idet ++){
                resp_arr[idet] += 1.0;
            }
            double sum = 0.0;
            for(int idet = 0; idet < ndet; idet ++){
                sum += resp_arr[idet];
            }
            for(int idet = 0; idet < ndet; idet ++){
                resp_arr[idet] /= sum;
            }

            int status = 0;
            int bitpix = -32;
            int naxis = 2;
            long* naxes = new long[naxis];
            naxes[0] = ndetx;
            naxes[1] = ndety;
            long npix_image = naxes[0] * naxes[1];
            fitsfile* fptr_out = NULL;
            fits_create_file(&fptr_out, outfile, &status);
            fits_create_img(fptr_out, bitpix, naxis, const_cast<long*>(naxes), &status);
            long firstpix[2] = {1,1};
            fits_write_pix(fptr_out, TDOUBLE, firstpix,
                           npix_image, const_cast<double*>(resp_arr), &status);
            fits_close_file(fptr_out, &status);
            fits_report_error(stderr, status);

            delete [] naxes;
            delete [] resp_arr;
        }
    }
    
    return status_prog;
}
