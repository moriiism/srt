#include "mi_sort.h"
#include "mir_math.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "mir_hist_info.h"
#include "mir_hist1d_nerr.h"
#include "mir_qdp_tool.h"
#include "arg_imgproj.h"
#include "sub.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    ArgValImgproj* argval = new ArgValImgproj;
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
    delete img_info_in;
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
    delete [] in_arr;


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
    sprintf(outqdp, "%s/%s_fine.qdp",
            argval->GetOutdir().c_str(),
            argval->GetOutfileHead().c_str());
    MirQdpTool::MkQdp(hd1d_proj_fine, outqdp, "x,y");
    
    for(int isky = 0; isky < nsky; isky++){
        // center
        Vect2d* vect = new Vect2d;
        vect->Init(hd2d->GetHi2d()->GetBinCenterXFromIbin(isky),
                   hd2d->GetHi2d()->GetBinCenterYFromIbin(isky));
        Vect2d* vect_conv = MirGeom::GenMotion(vect,
                                               argval->GetXpos(),
                                               argval->GetYpos(),
                                               theta_rad, -1);

        int ncor = 4;
        Vect2d** vect_c      = new Vect2d* [ncor];
        Vect2d** vect_c_conv = new Vect2d* [ncor];
        Vect2d** vect_c_conv_sort = new Vect2d* [ncor];
        for(int icor = 0; icor < ncor; icor++){
            vect_c[icor] = new Vect2d;
            vect_c_conv_sort[icor] = new Vect2d;
        }
        int icor = 0;
        vect_c[icor]->Init(hd2d->GetHi2d()->GetBinCenterXFromIbin(isky) - hd2d->GetHi2d()->GetBinWidthX() / 2.,
                           hd2d->GetHi2d()->GetBinCenterYFromIbin(isky) - hd2d->GetHi2d()->GetBinWidthY() / 2.);
        icor ++;
        vect_c[icor]->Init(hd2d->GetHi2d()->GetBinCenterXFromIbin(isky) + hd2d->GetHi2d()->GetBinWidthX() / 2.,
                           hd2d->GetHi2d()->GetBinCenterYFromIbin(isky) - hd2d->GetHi2d()->GetBinWidthY() / 2.);
        icor ++;
        vect_c[icor]->Init(hd2d->GetHi2d()->GetBinCenterXFromIbin(isky) + hd2d->GetHi2d()->GetBinWidthX() / 2.,
                           hd2d->GetHi2d()->GetBinCenterYFromIbin(isky) + hd2d->GetHi2d()->GetBinWidthY() / 2.);
        icor ++;
        vect_c[icor]->Init(hd2d->GetHi2d()->GetBinCenterXFromIbin(isky) - hd2d->GetHi2d()->GetBinWidthX() / 2.,
                           hd2d->GetHi2d()->GetBinCenterYFromIbin(isky) + hd2d->GetHi2d()->GetBinWidthY() / 2.);
        icor ++;
        for(int icor = 0; icor < ncor; icor++){
            vect_c_conv[icor] = MirGeom::GenMotion(vect_c[icor],
                                                   argval->GetXpos(),
                                                   argval->GetYpos(),
                                                   theta_rad, -1);
          
        }
        double* xpos_c_arr = new double [ncor];
        for(int icor = 0; icor < ncor; icor++){
            xpos_c_arr[icor] = vect_c_conv[icor]->GetPosX();
        }
        int index_xpos_c_arr[4];
        MiSort::Sort<double, int>(ncor, xpos_c_arr, index_xpos_c_arr, 0);
        delete [] xpos_c_arr;
        for(int icor = 0; icor < ncor; icor++){
            vect_c_conv_sort[icor]->Copy(vect_c_conv[ index_xpos_c_arr[icor] ]);
        }
        if(vect_c_conv_sort[0]->GetPosX() < argval->GetXloNew() ||
           vect_c_conv_sort[0]->GetPosX() > argval->GetXupNew() ||
           vect_c_conv_sort[1]->GetPosX() < argval->GetXloNew() ||
           vect_c_conv_sort[1]->GetPosX() > argval->GetXupNew() ||
           vect_c_conv_sort[2]->GetPosX() < argval->GetXloNew() ||
           vect_c_conv_sort[2]->GetPosX() > argval->GetXupNew() ||
           vect_c_conv_sort[3]->GetPosX() < argval->GetXloNew() ||
           vect_c_conv_sort[3]->GetPosX() > argval->GetXupNew()   ){
            continue;
        }
        double value_bin = hd2d->GetOvalElm(hd2d->GetHi2d()->GetIbinX(isky),
                                            hd2d->GetHi2d()->GetIbinY(isky));
        double area_bin  = hd2d->GetHi2d()->GetBinArea();
        FillTriLo(value_bin, area_bin,
                  vect_c_conv_sort[0],
                  vect_c_conv_sort[1],
                  vect_c_conv_sort[2],
                  argval->GetYLo(),
                  argval->GetYUp(),
                  hd1d_proj);
        FillTriUp(value_bin, area_bin,
                  vect_c_conv_sort[1],
                  vect_c_conv_sort[2],
                  vect_c_conv_sort[3],
                  argval->GetYLo(),
                  argval->GetYUp(),
                  hd1d_proj);
        FillSq(value_bin, area_bin,
               vect_c_conv_sort[0],
               vect_c_conv_sort[1],
               vect_c_conv_sort[2],
               vect_c_conv_sort[3],
               argval->GetYLo(),
               argval->GetYUp(),
               hd1d_proj);
        delete vect;
        delete vect_conv;
        for(int icor = 0; icor < ncor; icor++){
            delete vect_c[icor];
            delete vect_c_conv[icor];
            delete vect_c_conv_sort[icor];
        }
        delete [] vect_c;
        delete [] vect_c_conv;
        delete [] vect_c_conv_sort;
    }

    MirQdpTool::MkQdp(hd1d_proj, "temp.qdp", "x,y");

    delete hd2d;
    delete hd1d_proj;
    delete argval;
    
    return status_prog;
}


