#include "mir_math.h"
#include "mib_blas.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "arg_calc_flux0.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    ArgValCalcFlux0* argval = new ArgValCalcFlux0;
    argval->Init(argc, argv);
    argval->Print(stdout);

    // load data_on_file
    int ndetx = MifFits::GetAxisSize(argval->GetDataOnFile(), 0);
    int ndety = MifFits::GetAxisSize(argval->GetDataOnFile(), 1);
    int ndet = ndetx * ndety;
    double* data_on_arr = NULL;
    MifImgInfo* img_info_data_on = new MifImgInfo;
    img_info_data_on->InitSetImg(1, 1, ndetx, ndety);
    int bitpix_data_on = 0;
    MifFits::InFitsImageD(argval->GetDataOnFile(),
                          img_info_data_on,
                          &bitpix_data_on, &data_on_arr);
    double nevt_on = MirMath::GetSum(ndet, data_on_arr);
    printf("nevt_on = %e\n", nevt_on);
    delete [] data_on_arr;
    delete img_info_data_on;

    // load data_off_file
    double* data_off_arr = NULL;
    MifImgInfo* img_info_data_off = new MifImgInfo;
    img_info_data_off->InitSetImg(1, 1, ndetx, ndety);
    int bitpix_data_off = 0;
    MifFits::InFitsImageD(argval->GetDataOffFile(),
                          img_info_data_off,
                          &bitpix_data_off, &data_off_arr);
    double nevt_off = MirMath::GetSum(ndet, data_off_arr);
    printf("nevt_off = %e\n", nevt_off);
    delete [] data_off_arr;
    delete img_info_data_off;

    double phase_ratio_on = argval->GetPhaseRatioOn();
    double phase_ratio_off = argval->GetPhaseRatioOff();
    double live_time_ratio_on = argval->GetLiveTimeRatioOn();
    double live_time_ratio_off = argval->GetLiveTimeRatioOff();

    double flux0 = nevt_on /(phase_ratio_on * live_time_ratio_on)
        - nevt_off /(phase_ratio_off * live_time_ratio_off);
    printf("flux0 = %e\n", flux0);
    
    return status_prog;
}
