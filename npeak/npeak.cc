#include "mi_sort.h"
#include "mir_math.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "mir_hist_info.h"
#include "mir_hist1d_nerr.h"
#include "mir_qdp_tool.h"
#include "arg_npeak.h"
#include "sub.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    ArgValNpeak* argval = new ArgValNpeak;
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
    int nsky = nskyx * nskyy;
    
    // load infile
    double* in_arr = NULL;
    MifImgInfo* img_info_in = new MifImgInfo;
    img_info_in->InitSetImg(1, 1, nskyx, nskyy);
    int bitpix_in = 0;
    MifFits::InFitsImageD(argval->GetInfile(), img_info_in,
                          &bitpix_in, &in_arr);
    double nph_in = MirMath::GetSum(nsky, in_arr);
    printf("nph_in = %e\n", nph_in);

    HistDataNerr1d* hd1d_proj = new HistDataNerr1d;
    hd1d_proj->Init(argval->GetNbinxNew(), argval->GetXloNew(), argval->GetXupNew());
    HistDataNerr1d* hd1d_proj_fine = new HistDataNerr1d;
    hd1d_proj_fine->Init(argval->GetNbinxNew(), argval->GetXloNew(), argval->GetXupNew());
    double theta_rad = argval->GetTheta() / 180. * M_PI;
    
    // 2d hist of input image
    HistDataNerr2d* hd2d = new HistDataNerr2d;
    hd2d->Init(nskyx, 0.0, nskyx,
               nskyy, 0.0, nskyy);
    hd2d->SetOvalArr(nsky, in_arr);

    // 2d hist of input image with fine resolution
    int factor_fine = 10;
    int nskyx_fine = nskyx * factor_fine;
    int nskyy_fine = nskyy * factor_fine;
    int nsky_fine = nskyx_fine * nskyy_fine;
    HistDataNerr2d* hd2d_fine = new HistDataNerr2d;
    hd2d_fine->Init(nskyx_fine, 0.0, nskyx,
                    nskyy_fine, 0.0, nskyy);
    for(int isky_fine = 0; isky_fine < nsky_fine; isky_fine ++){
        double posx = hd2d_fine->GetHi2d()->GetBinCenterXFromIbin(isky_fine);
        double posy = hd2d_fine->GetHi2d()->GetBinCenterYFromIbin(isky_fine);
        int isky = hd2d->GetHi2d()->GetIbinFromXY(posx, posy);
        hd2d_fine->Fill(posx, posy, in_arr[isky] / pow(factor_fine, 2) );
    }

 

    // fine
    for(int isky_fine = 0; isky_fine < nsky_fine; isky_fine ++){
        // center
        Vect2d* vect = new Vect2d;
        vect->Init(hd2d_fine->GetHi2d()->GetBinCenterXFromIbin(isky_fine),
                   hd2d_fine->GetHi2d()->GetBinCenterYFromIbin(isky_fine));
        Vect2d* vect_conv = MirGeom::GenMotion(vect,
                                               argval->GetXpos(),
                                               argval->GetYpos(),
                                               theta_rad, -1);
        if(vect_conv->GetPosX() < argval->GetXloNew() ||
           vect_conv->GetPosX() > argval->GetXupNew()   ){
            continue;
        }
        if(argval->GetYLo() <= vect_conv->GetPosY() &&
           vect_conv->GetPosY() <= argval->GetYUp()   ){
            hd1d_proj_fine->Fill(vect_conv->GetPosX(),
                                 hd2d_fine->GetOvalArr()->GetValElm(isky_fine));
        }

    }
    char outqdp[kLineSize];
    sprintf(outqdp, "%s/%s_proj.qdp",
            argval->GetOutdir().c_str(),
            argval->GetOutfileHead().c_str());
    MirQdpTool::MkQdp(hd1d_proj_fine, outqdp, "x,y");

    //
    // calc noise level
    //
    HistDataNerr2d* hd2d_mask_bg = new HistDataNerr2d;
    hd2d_mask_bg->Copy(hd2d);
    hd2d_mask_bg->SetConst(1.0);
    
    double mean   = 0.0;
    double stddev = 0.0;
    int nclip = 30;
    double significance_clip = 3.0;
    for(int iclip = 0; iclip < nclip; iclip++){
        // get mean, stddev
        vector<double> val_vec;
        for(long ibin = 0; ibin < hd2d->GetNbin(); ibin++){
            val_vec.push_back(hd2d->GetOvalArr()->GetValElm(ibin));
        }
        mean   = MirMath::GetAMean(val_vec);
        stddev = MirMath::GetSqrtOfUnbiasedVariance(val_vec);
        //        printf("iclip, mean, stddev, nsize = %d, %e, %e, %d\n",
        //               iclip, mean, stddev, (int) val_vec.size());

        for(long ibin = 0; ibin < hd2d->GetNbin(); ibin++){
            if( hd2d_mask_bg->GetOvalArr()->GetValElm(ibin) > 0){
                if( fabs(hd2d->GetOvalArr()->GetValElm(ibin) - mean) >  significance_clip * stddev){
                    hd2d_mask_bg->SetOvalElm(hd2d_mask_bg->GetHi2d()->GetIbinX(ibin),
                                             hd2d_mask_bg->GetHi2d()->GetIbinY(ibin),
                                             0);
                } else{
                    hd2d_mask_bg->SetOvalElm(hd2d_mask_bg->GetHi2d()->GetIbinX(ibin),
                                             hd2d_mask_bg->GetHi2d()->GetIbinY(ibin),
                                             1);
                }
            }
        }

        long nbin_mask_bg = 0;
        for(long ibin = 0; ibin < hd2d->GetNbin(); ibin++){
            if( hd2d_mask_bg->GetOvalArr()->GetValElm(ibin) > 0){
                nbin_mask_bg ++;
            }
        }
        if(nbin_mask_bg == (int) val_vec.size()){
            break;
        }
    }

    
//    MifFits::OutFitsImageD(argval->GetOutdir(),
//                           argval->GetOutfileHead(),
//                           "mask_bg", 2,
//                           bitpix_in,
//                           img_info_in->GetNaxesArr(),
//                           hd2d_mask_bg->GetOvalArr()->GetVal());

    
    double width = argval->GetYUp() - argval->GetYLo();
    printf("width = %e\n", width);
    printf("mean, stddev = %e, %e\n", mean, stddev);
    printf("mean * width, stddev * width = %e, %e\n", mean * width, stddev * width);
    
    double threshold = (mean + argval->GetSignificance() * stddev) * width;
    printf("threshold = %e\n", threshold);
    
    Interval* interval = new Interval;
    vector<double> tstart_vec;
    vector<double> tstop_vec;
    for(int ibin = 0; ibin < hd1d_proj_fine->GetNbinX(); ibin ++){
        if(hd1d_proj_fine->GetOvalElm(ibin) > threshold){
            tstart_vec.push_back(hd1d_proj_fine->GetBinLo(ibin));
            tstop_vec.push_back(hd1d_proj_fine->GetBinUp(ibin));
        }
    }
    interval->Init(tstart_vec.size());
    interval->Set(tstart_vec, tstop_vec);
    interval->Clean(hd1d_proj_fine->GetXvalBinWidth() / 10.0);

    char outqdp_src[kLineSize];
    sprintf(outqdp_src, "%s/%s_proj_src.qdp",
            argval->GetOutdir().c_str(),
            argval->GetOutfileHead().c_str());
    if( 0 < tstart_vec.size() ){
        MirQdpTool::MkQdp(interval, outqdp_src);
    }

    char outdat[kLineSize];
    sprintf(outdat, "%s/%s_out.dat",
            argval->GetOutdir().c_str(),
            argval->GetOutfileHead().c_str());
    FILE* fp_outdat = fopen(outdat, "w");
    fprintf(fp_outdat, "%ld\n", interval->GetNterm());
    printf("%ld\n", interval->GetNterm());
    fclose(fp_outdat);

    
    
    delete hd2d;
    delete hd1d_proj;
    delete argval;
    
    return status_prog;
}


