#include "mir_math.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "arg_mkbg.h"
#include "model.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    ArgValMkbg* argval = new ArgValMkbg;
    argval->Init(argc, argv);
    argval->Print(stdout);

    int nskyx = 60;
    int nskyy = 60;
    int nsky = nskyx * nskyy;
    double* sky_arr = new double[nsky];
    for(int isky = 0; isky < nsky; isky ++){
        sky_arr[isky] = 0.0;
    }
    
    // load params of model functions
    string* lines_arr = NULL;
    long nline = 0;
    MiIolib::GenReadFileSkipComment(argval->GetModelFile(),
                                    &lines_arr, &nline);
    if (nline != 2){
        abort();
    }
    string model  = lines_arr[0];
    string params = lines_arr[1];
    string* par_arr = NULL;
    int npar = 0;
    MiStr::GenSplit(params, &npar, &par_arr);
    
    if (model == "const"){
        if(npar != 1){
            abort();
        }
        double constant = atof(par_arr[0].c_str());
        printf("%e\n", constant);
        double val_min = 10.0;
        for(int isky = 0; isky < nsky; isky ++){
            int posx = isky % nskyx;
            int posy = isky / nskyx;
            double val = GetConstFunc(posx, posy, constant);
            sky_arr[isky] = val;
            if (val < val_min){
                val_min = val;
            }
        }
        if(val_min < 0.0){
            printf("Error: val_min(%e) < 0.0, then abort().\n", val_min);
            printf("Add val_min on constant.\n");
            abort();
        }
    } else if(model == "lin"){
        if(npar != 3){
            abort();
        }
        double w0 = atof(par_arr[0].c_str());
        double w1 = atof(par_arr[1].c_str());
        double constant = atof(par_arr[2].c_str());
        printf("%e  %e  %e\n",
               w0, w1, constant);
        double val_min = 10.0;
        for(int isky = 0; isky < nsky; isky ++){
            int posx = isky % nskyx;
            int posy = isky / nskyx;
            double val = GetLinFunc(posx, posy,
                                    w0, w1, constant);
            sky_arr[isky] = val;
            if (val < val_min){
                val_min = val;
            }
        }
        if(val_min < 0.0){
            printf("Error: val_min(%e) < 0.0, then abort().\n", val_min);
            printf("Add val_min on constant.\n");
            abort();
        }
    } else if(model == "quad"){
        if(npar != 6){
            abort();
        }
        double posx_ctr = atof(par_arr[0].c_str());
        double posy_ctr = atof(par_arr[1].c_str());
        double coeff0 = atof(par_arr[2].c_str());
        double coeff1 = atof(par_arr[3].c_str());
        double coeff2 = atof(par_arr[4].c_str());
        double constant = atof(par_arr[5].c_str());
        printf("%e  %e  %e  %e  %e  %e\n",
               posx_ctr, posy_ctr, coeff0, coeff1, coeff2, constant);
        double val_min = 10.0;
        for(int isky = 0; isky < nsky; isky ++){
            int posx = isky % nskyx;
            int posy = isky / nskyx;
            double val = GetQuadFunc(posx, posy,
                                     posx_ctr, posy_ctr,
                                     coeff0, coeff1, coeff2,
                                     constant);
            sky_arr[isky] = val;
            if (val < val_min){
                val_min = val;
            }
        }
        if(val_min < 0.0){
            printf("Error: val_min(%e) < 0.0, then abort().\n", val_min);
            printf("Add val_min on constant.\n");
            abort();
        }
    } else{
        abort();
    }
    delete [] par_arr;
    
    // normalize
    double sum = 0.0;
    for(int isky = 0; isky < nsky; isky ++){
        sum += sky_arr[isky];
    }
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
