#include "mir_math.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "arg_mkimg.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    ArgValMkimg* argval = new ArgValMkimg;
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
    
    // load infile
    double* sky_arr = NULL;
    MifImgInfo* img_info_in = new MifImgInfo;
    img_info_in->InitSetImg(1, 1, nskyx, nskyy);
    int bitpix_in = 0;
    MifFits::InFitsImageD(argval->GetInfile(), img_info_in,
                          &bitpix_in, &sky_arr);

    string* lines_arr = NULL;
    long nline = 0;
    MiIolib::GenReadFileSkipComment(argval->GetDotfile(),
                                    &lines_arr, &nline);
    long* xpos_arr = new long [nline];
    long* ypos_arr = new long [nline];
    double* flux_arr = new double [nline];
    for(long iline = 0; iline < nline; iline ++){
        int nsplit = 0;
        string* split_arr = NULL;
        MiStr::GenSplit(lines_arr[iline], &nsplit, &split_arr);
        xpos_arr[iline] = atoi(split_arr[0].c_str());
        ypos_arr[iline] = atoi(split_arr[1].c_str());
        flux_arr[iline] = atof(split_arr[2].c_str());
        delete [] split_arr;
    }
    MiIolib::DelReadFile(lines_arr);

    for(long index = 0; index < nskyx * nskyy; index ++){
        sky_arr[index] = 0.0;
    }
    
    for(long iline = 0; iline < nline; iline ++){
        long index = ypos_arr[iline] * img_info_in->GetNaxesArrElm(0) + xpos_arr[iline];
        sky_arr[index] = flux_arr[iline];
    }
    
    {
        int naxis = 2;
        long* naxes = new long[naxis];
        naxes[0] = nskyx;
        naxes[1] = nskyy;
        char tag[kLineSize];
        sprintf(tag, "mkimg");
        int bitpix = -64;
        MifFits::OutFitsImageD(argval->GetOutdir(), argval->GetOutfileHead(),
                               tag, 2,
                               bitpix,
                               naxes, sky_arr);
        delete [] naxes;
    }
    


    return status_prog;
}
