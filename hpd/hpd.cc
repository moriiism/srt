#include "mir_math.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "mir_hist_info.h"
#include "mir_hist1d_nerr.h"
#include "mir_qdp_tool.h"
#include "arg_hpd.h"
#include "sub.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    ArgValHpd* argval = new ArgValHpd;
    argval->Init(argc, argv);
    argval->Print(stdout);


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
    delete img_info_in;
    double nph_in = MirMath::GetSum(nsky, in_arr);
    printf("nph_in = %e\n", nph_in);

    // 2d hist_info of input image
    HistInfo2d* hi2d = new HistInfo2d;
    hi2d->InitSetByNbin(0.0, nskyx, nskyx, 0.0, nskyy, nskyy);

    double xpos0 = argval->GetXpos();
    double ypos0 = argval->GetYpos();
    HistDataNerr1d* hd1d = new HistDataNerr1d;
    hd1d->Init(100, 0.0, 60);
    for(int isky = 0; isky < nsky; isky++){
        double xpos = hi2d->GetBinCenterXFromIbin(isky);
        double ypos = hi2d->GetBinCenterYFromIbin(isky);
        double rad = sqrt( pow(xpos - xpos0, 2) + pow(ypos - ypos0, 2) );
        hd1d->Fill(rad, in_arr[isky]);
    }

    HistDataNerr1d* hd1d_cum = new HistDataNerr1d;
    hd1d_cum->Init(100, 0.0, 60);
    int nrad = 100;
    double sum = 0.0;
    for(int irad = 0; irad < nrad; irad++){
        sum += hd1d->GetOvalElm(irad);
        hd1d_cum->SetOvalElm(irad, sum);
    }
    MirQdpTool::MkQdp(hd1d, "temp.qdp", "x,y");
    MirQdpTool::MkQdp(hd1d_cum, "sum.qdp", "x,y");

    double half_val = hd1d_cum->GetOvalElm(nrad - 1) / 2.;
    double hpr = 0.0;
    for(int irad = 0; irad < nrad; irad++){
        if(half_val < hd1d_cum->GetOvalElm(irad)){
            printf("irad = %d\n", irad);
            hpr = irad;
            break;
        }
    }
    double hpd = hpr * 2.0;
    printf("hpd = %e\n", hpd);
    
    return status_prog;
}
