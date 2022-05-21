#include "mir_math.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "arg_mksrc.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    ArgValMksrc* argval = new ArgValMksrc;
    argval->Init(argc, argv);
    argval->Print(stdout);

   
    string* lines_arr = NULL;
    long nline = 0;
    MiIolib::GenReadFileSkipComment(argval->GetInfile(), &lines_arr, &nline);
    int npos = nline;
    int* xpos_arr = new int [npos];
    int* ypos_arr = new int [npos];
    double* val_arr  = new double [npos];

    for(long iline = 0; iline < nline; iline ++){
        string* split_arr = NULL;
        int nsplit = 0;
        MiStr::GenSplit(lines_arr[iline], &nsplit, &split_arr);
        xpos_arr[iline] = atoi(split_arr[0].c_str());
        ypos_arr[iline] = atoi(split_arr[1].c_str());
        val_arr[iline]  = atof(split_arr[2].c_str());
        printf("%d  %d  %e\n", xpos_arr[iline], ypos_arr[iline], val_arr[iline]);
        delete [] split_arr;
        
    }
    
    int nskyx = argval->GetNskyx();
    int nskyy = argval->GetNskyy();
    int nsky = nskyx * nskyy;
    double* sky_arr = new double[nsky];
    for(int isky = 0; isky < nsky; isky ++){
        sky_arr[isky] = 0.0;
    }

    double sum = 0.0;
    for(int ipos = 0; ipos < npos; ipos ++){
        int isky_pos = ypos_arr[ipos] * nskyx + xpos_arr[ipos];
        sky_arr[isky_pos] = val_arr[ipos];
        sum += val_arr[ipos];
    }
    // normalize
    for(int isky = 0; isky < nsky; isky ++){
        sky_arr[isky] /= sum;
    }

    int status = 0;
    int naxis = 2;
    long* naxes = new long[naxis];
    naxes[0] = nskyx;
    naxes[1] = nskyy;    
    long npix_image = naxes[0] * naxes[1];
    int bitpix = -64;
    fitsfile* fptr_out = NULL;
    fits_create_file(&fptr_out, argval->GetOutfile().c_str(), &status);
    fits_create_img(fptr_out, bitpix, naxis, const_cast<long*>(naxes), &status);

    long firstpix[2] = {1,1};
    fits_write_pix(fptr_out, TDOUBLE, firstpix,
                   npix_image, const_cast<double*>(sky_arr), &status);
    fits_close_file(fptr_out, &status);
    fits_report_error(stderr, status);
    

    
    return status_prog;
}
