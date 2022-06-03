#include "mir_math.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "arg_mkimg_points.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    ArgValMkimgPoints* argval = new ArgValMkimgPoints;
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
    double* sky_norm_arr = new double[nsky];
    for(int isky = 0; isky < nsky; isky ++){
        sky_arr[isky] = 0.0;
        sky_norm_arr[isky] = 0.0;
    }

    double sum = 0.0;
    for(int ipos = 0; ipos < npos; ipos ++){
        int isky_pos = ypos_arr[ipos] * nskyx + xpos_arr[ipos];
        sky_arr[isky_pos] = val_arr[ipos];
        sum += val_arr[ipos];
    }
    // normalize
    for(int isky = 0; isky < nsky; isky ++){
        sky_norm_arr[isky] = sky_arr[isky] / sum;
    }

    long naxes[2];
    naxes[0] = nskyx;
    naxes[1] = nskyy;
    int bitpix = -32;
    char outfile[kLineSize];
    sprintf(outfile, "%s/%s.fits", argval->GetOutdir().c_str(),
            argval->GetOutfileHead().c_str());
    MifFits::OutFitsImageD(outfile, 2, bitpix,
                           naxes, sky_arr);

    sprintf(outfile, "%s/%s_norm.fits", argval->GetOutdir().c_str(),
            argval->GetOutfileHead().c_str());
    MifFits::OutFitsImageD(outfile, 2, bitpix,
                           naxes, sky_norm_arr);
    
    return status_prog;
}
