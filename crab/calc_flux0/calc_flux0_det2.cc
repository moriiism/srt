#include "mib_blas.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "srtmathlib.h"
#include "arg_calc_flux0_det2.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    ArgValCalcFlux0Det2* argval = new ArgValCalcFlux0Det2;
    argval->Init(argc, argv);
    argval->Print(stdout);

    // load data_on_file
    int ndetx = MifFits::GetAxisSize(argval->GetDataOnFileDet1(), 0);
    int ndety = MifFits::GetAxisSize(argval->GetDataOnFileDet1(), 1);
    int ndet = ndetx * ndety;

    double* data_on_det1_arr = NULL;
    double* data_on_det2_arr = NULL;
    MifImgInfo* img_info_data_on = new MifImgInfo;
    img_info_data_on->InitSetImg(1, 1, ndetx, ndety);
    int bitpix_data_on = 0;
    MifFits::InFitsImageD(argval->GetDataOnFileDet1(),
                          img_info_data_on,
                          &bitpix_data_on, &data_on_det1_arr);
    MifFits::InFitsImageD(argval->GetDataOnFileDet2(),
                          img_info_data_on,
                          &bitpix_data_on, &data_on_det2_arr);
    double nevt_on_det1 = SrtMathlib::GetSum(ndet, data_on_det1_arr);
    double nevt_on_det2 = SrtMathlib::GetSum(ndet, data_on_det2_arr);
    printf("nevt_on_det1 = %e\n", nevt_on_det1);
    printf("nevt_on_det2 = %e\n", nevt_on_det2);
    delete [] data_on_det1_arr;
    delete [] data_on_det2_arr;
    delete img_info_data_on;

    
    double nevt_off_det1 = 0.0;
    double nevt_off_det2 = 0.0;
    if(argval->GetDataOffFileDet1() != "none"){
        // load data_off_file
        double* data_off_det1_arr = NULL;
        double* data_off_det2_arr = NULL;
        MifImgInfo* img_info_data_off = new MifImgInfo;
        img_info_data_off->InitSetImg(1, 1, ndetx, ndety);
        int bitpix_data_off = 0;
        MifFits::InFitsImageD(argval->GetDataOffFileDet1(),
                              img_info_data_off,
                              &bitpix_data_off, &data_off_det1_arr);
        MifFits::InFitsImageD(argval->GetDataOffFileDet2(),
                              img_info_data_off,
                              &bitpix_data_off, &data_off_det2_arr);
        nevt_off_det1 = SrtMathlib::GetSum(ndet, data_off_det1_arr);
        nevt_off_det2 = SrtMathlib::GetSum(ndet, data_off_det2_arr);
        printf("nevt_off_det1 = %e\n", nevt_off_det1);
        printf("nevt_off_det2 = %e\n", nevt_off_det2);
        delete [] data_off_det1_arr;
        delete [] data_off_det2_arr;
        delete img_info_data_off;
    } else {
        printf("data_off_file == none, then flux_off = 0\n");
    }

    double live_time_ratio_on_det1
        = argval->GetLiveTimeRatioOnDet1();
    double live_time_ratio_on_det2
        = argval->GetLiveTimeRatioOnDet2();
    double live_time_ratio_off_det1
        = argval->GetLiveTimeRatioOffDet1();
    double live_time_ratio_off_det2
        = argval->GetLiveTimeRatioOffDet2();    
    double phase_ratio_on = argval->GetPhaseRatioOn();
    double phase_ratio_off = argval->GetPhaseRatioOff();

    double flux_on = (
        nevt_on_det1 / live_time_ratio_on_det1 +
        nevt_on_det2 / live_time_ratio_on_det2) /
        phase_ratio_on;
    
    double flux0 = 0.0;
    double flux0_err = 0.0;
    if(argval->GetDataOffFileDet1() != "none"){
        double flux_off = (
            nevt_off_det1 / live_time_ratio_off_det1 +
            nevt_off_det2 / live_time_ratio_off_det2) /
            phase_ratio_off;
        flux0 = flux_on - flux_off;
        flux0_err = sqrt(
            nevt_on_det1 / pow(
                phase_ratio_on *
                live_time_ratio_on_det1, 2) +
            nevt_on_det2 / pow(
                phase_ratio_on *
                live_time_ratio_on_det2, 2) +
            nevt_off_det1 / pow(
                phase_ratio_off *
                live_time_ratio_off_det1, 2) +
            nevt_off_det2 / pow(
                phase_ratio_off *
                live_time_ratio_off_det2, 2));
    } else {
        flux0 = flux_on;
        flux0_err = sqrt(
            nevt_on_det1 / pow(
                phase_ratio_on *
                live_time_ratio_on_det1, 2) +
            nevt_on_det2 / pow(
                phase_ratio_on *
                live_time_ratio_on_det2, 2));
    }
    printf("flux0  flux0_err\n");
    printf("%e  %e\n", flux0, flux0_err);

    return status_prog;
}
