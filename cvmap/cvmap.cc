#include "mi_iolib.h"
#include "mir_root_tool.h"
#include "mir_hist2d_nerr.h"
#include "mir_qdp_tool.h"
#include "arg_cvmap.h"
#include "TColor.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[]){
    int status = kRetNormal;
  
    ArgValCvmap* argval = new ArgValCvmap;
    argval->Init(argc, argv);

    MirRootTool* root_tool = new MirRootTool;
    root_tool->InitTCanvas("pub");
    if( MiIolib::TestFileExist(argval->GetOutdir()) ){
        char cmd[kLineSize];
        sprintf(cmd, "mkdir -p %s", argval->GetOutdir().c_str());
        system(cmd);
    }


    string* lines_arr = NULL;
    long nline = 0;
    MiIolib::GenReadFileSkipComment(argval->GetInfile(),
                                    &lines_arr,
                                    &nline);
    double* imu_arr = new double [nline];
    double* ibeta_arr = new double [nline];
    double* mu_arr = new double [nline];
    double* beta_arr = new double [nline];
    double* ave_arr = new double [nline];
    double* stddev_arr = new double [nline];
    for(long iline = 0; iline < nline; iline ++){
        int nsplit = 0;
        string* split_arr = NULL;
        MiStr::GenSplit(lines_arr[iline], &nsplit, &split_arr);

        imu_arr[iline]   = atof(split_arr[0].c_str());
        ibeta_arr[iline] = atof(split_arr[1].c_str());
        mu_arr[iline]    = atof(split_arr[2].c_str());
        beta_arr[iline]  = atof(split_arr[3].c_str());
        ave_arr[iline]   = atof(split_arr[4].c_str());
        stddev_arr[iline] = atof(split_arr[5].c_str());
        MiStr::DelSplit(split_arr);
        printf("%e %e  %e %e  %e %e\n",
               imu_arr[iline], ibeta_arr[iline],
               mu_arr[iline], beta_arr[iline],
               ave_arr[iline], stddev_arr[iline]);
    }

    double ave_min   = 1.0e10;
    double ave_max   = -1.0e10;
    double imu_min   = 0.0;
    double ibeta_min = 0.0;
    double mu_min    = 0.0;
    double beta_min  = 0.0;
    HistInfo2d* hi2d = new HistInfo2d;
    hi2d->Load(argval->GetHistInfoFile());
    HistDataNerr2d* hd2d = new HistDataNerr2d;
    hd2d->Init(hi2d);
    HistDataNerr2d* hd2d_stddev = new HistDataNerr2d;
    hd2d_stddev->Init(hi2d);
    for(long iline = 0; iline < nline; iline ++){
        hd2d->Fill(imu_arr[iline], ibeta_arr[iline], ave_arr[iline]);
        hd2d_stddev->Fill(imu_arr[iline], ibeta_arr[iline], stddev_arr[iline]);

        if(ave_arr[iline] < ave_min){
            ave_min   = ave_arr[iline];
            imu_min   = imu_arr[iline];
            ibeta_min = ibeta_arr[iline];
            mu_min    = mu_arr[iline];
            beta_min  = beta_arr[iline];            
        }
        if(ave_arr[iline] > ave_max){
            ave_max   = ave_arr[iline];
        }
        
    }

    //printf("min at (%e, %e)\n",
    //       hd2d->GetXvalAtOvalMin(),
    //       hd2d->GetYvalAtOvalMin());
    printf("min at (imu, ibeta) (mu, beta) = (%e, %e) (%e, %e)\n",
           imu_min, ibeta_min, mu_min, beta_min);

    printf("ave_min = %e, ave_max = %e\n", ave_min, ave_max);
    
    double zrange_lo = argval->GetZrangeLo();
    double zrange_up = argval->GetZrangeUp();
    
    TH2D* th2d = hd2d->GenTH2D(0, 0, 0);
    th2d->SetAxisRange(zrange_lo, zrange_up, "Z");
    gStyle->SetPalette(53);
    // TColor::InvertPalette();
    th2d->Draw("COLZ");
    
    
    gPad->Update();
    // TPaletteAxis* palette = (TPaletteAxis*) th2d->GetListOfFunctions()->FindObject("palette");
//    palette->SetX1NDC(0.86);
//    palette->SetX2NDC(0.89);
    th2d->GetXaxis()->SetTitleSize(0.05);
    th2d->GetYaxis()->SetTitleSize(0.05);
    th2d->GetXaxis()->SetLabelSize(0.05);
    th2d->GetYaxis()->SetLabelSize(0.05);

    char outfig[kLineSize];
    sprintf(outfig, "%s/%s.png",
            argval->GetOutdir().c_str(),
            argval->GetOutfileHead().c_str());
//    hd2d->MkTH2FigZrange(outfig, root_tool,
//                         2.27e-2, 2.300e-2,
//                         //7.2e-2, 7.3e-2,
//                         0.0, 0.0, 0.0,
//                         "mu", "beta", "ave");

    root_tool->GetTCanvas()->Print(outfig);

    /// stddev
    TH2D* th2d_stddev = hd2d_stddev->GenTH2D(0, 0, 0);
    th2d_stddev->Draw("COLZ");
    gPad->Update();
    th2d_stddev->GetXaxis()->SetTitleSize(0.05);
    th2d_stddev->GetYaxis()->SetTitleSize(0.05);
    th2d_stddev->GetXaxis()->SetLabelSize(0.05);
    th2d_stddev->GetYaxis()->SetLabelSize(0.05);

    char outfig_stddev[kLineSize];
    sprintf(outfig_stddev, "%s/%s_stddev.png",
            argval->GetOutdir().c_str(),
            argval->GetOutfileHead().c_str());
    
    root_tool->GetTCanvas()->Print(outfig_stddev);

    delete argval;
    
    return status;
}


