#include "mir_math.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "arg_mkimg_diffuse.h"
#include "model.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    ArgValMkimgDiffuse* argval = new ArgValMkimgDiffuse;
    argval->Init(argc, argv);
    argval->Print(stdout);
   
    int nskyx = argval->GetNskyx();
    int nskyy = argval->GetNskyy();
    int nsky = nskyx * nskyy;
    double* sky_arr = new double[nsky];
    double* sky_norm_arr = new double[nsky];    
    for(int isky = 0; isky < nsky; isky ++){
        sky_arr[isky] = 0.0;
        sky_norm_arr[isky] = 0.0;        
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
        for(int isky = 0; isky < nsky; isky ++){
            int posx = isky % nskyx;
            int posy = isky / nskyx;
            double val = GetConstFunc(posx, posy, constant);
            sky_arr[isky] = MirMath::GetMax(val, 0.0);
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
        for(int isky = 0; isky < nsky; isky ++){
            int posx = isky % nskyx;
            int posy = isky / nskyx;
            double val = GetLinFunc(posx, posy,
                                    w0, w1, constant);
            sky_arr[isky] = MirMath::GetMax(val, 0.0);
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
        for(int isky = 0; isky < nsky; isky ++){
            int posx = isky % nskyx;
            int posy = isky / nskyx;
            double val = GetQuadFunc(posx, posy,
                                     posx_ctr, posy_ctr,
                                     coeff0, coeff1, coeff2,
                                     constant);
            sky_arr[isky] = MirMath::GetMax(val, 0.0);            
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
