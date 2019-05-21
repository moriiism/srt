#include "mi_iolib.h"
#include "mir_root_tool.h"
#include "mir_hist2d_nerr.h"
#include "mir_qdp_tool.h"
#include "arg_ptmap.h"
#include "TColor.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[]){
    int status = kRetNormal;
  
    ArgValPtmap* argval = new ArgValPtmap;
    argval->Init(argc, argv);

    MirRootTool* root_tool = new MirRootTool;
    root_tool->InitTCanvas("pub");
    if( MiIolib::TestFileExist(argval->GetOutdir()) ){
        char cmd[kLineSize];
        sprintf(cmd, "mkdir -p %s", argval->GetOutdir().c_str());
        system(cmd);
    }


    string* lines_list_arr = NULL;
    long nline_list = 0;
    MiIolib::GenReadFileSkipComment(argval->GetCvsmrList(),
                                    &lines_list_arr,
                                    &nline_list);
    double** imu_arr    = new double* [nline_list];
    double** ibeta_arr  = new double* [nline_list];
    double** mu_arr     = new double* [nline_list];
    double** beta_arr   = new double* [nline_list];
    double** ave_arr    = new double* [nline_list];
    double** stddev_arr = new double* [nline_list];
    int*     nline_arr  = new int [nline_list];
    
    for(int iline_list = 0; iline_list < nline_list; iline_list ++){
        string* lines_arr = NULL;
        long nline = 0;
        MiIolib::GenReadFileSkipComment(lines_list_arr[iline_list],
                                        &lines_arr,
                                        &nline);
        nline_arr[iline_list] = nline;
        imu_arr[iline_list]    = new double [nline];
        ibeta_arr[iline_list]  = new double [nline];
        mu_arr[iline_list]     = new double [nline];
        beta_arr[iline_list]   = new double [nline];
        ave_arr[iline_list]    = new double [nline];
        stddev_arr[iline_list] = new double [nline];
        
        for(long iline = 0; iline < nline; iline ++){
            int nsplit = 0;
            string* split_arr = NULL;
            MiStr::GenSplit(lines_arr[iline], &nsplit, &split_arr);

            imu_arr[iline_list][iline]    = atof(split_arr[0].c_str());
            ibeta_arr[iline_list][iline]  = atof(split_arr[1].c_str());
            mu_arr[iline_list][iline]     = atof(split_arr[2].c_str());
            beta_arr[iline_list][iline]   = atof(split_arr[3].c_str());
            ave_arr[iline_list][iline]    = atof(split_arr[4].c_str());
            stddev_arr[iline_list][iline] = atof(split_arr[5].c_str());
            MiStr::DelSplit(split_arr);
            printf("%e %e  %e %e  %e %e\n",
                   imu_arr[iline_list][iline], ibeta_arr[iline_list][iline],
                   mu_arr[iline_list][iline], beta_arr[iline_list][iline],
                   ave_arr[iline_list][iline], stddev_arr[iline_list][iline]);
        }
        delete [] lines_arr;
    }


    HistInfo2d* hi2d = new HistInfo2d;
    hi2d->Load(argval->GetHistInfoFile());

    HistData2d** hd2d_arr = new HistData2d* [nline_list];
    HistData2d** hd2d_stddev_arr = new HistData2d* [nline_list];
    for(int iline_list = 0; iline_list < nline_list; iline_list ++){
        hd2d_arr[iline_list] = new HistDataNerr2d;
        hd2d_arr[iline_list]->Init(hi2d);
        hd2d_stddev_arr[iline_list] = new HistDataNerr2d;
        hd2d_stddev_arr[iline_list]->Init(hi2d);
    }
    
    HistInfo2d* hi2d_index = new HistInfo2d;
    hi2d_index->Load(argval->GetHistInfoIndexFile());
    HistData2d** hd2d_index_arr = new HistData2d* [nline_list];
    HistData2d** hd2d_stddev_index_arr = new HistData2d* [nline_list];
    for(int iline_list = 0; iline_list < nline_list; iline_list ++){
        hd2d_index_arr[iline_list] = new HistDataNerr2d;
        hd2d_index_arr[iline_list]->Init(hi2d_index);
        hd2d_stddev_index_arr[iline_list] = new HistDataNerr2d;
        hd2d_stddev_index_arr[iline_list]->Init(hi2d_index);
    }
    for(int iline_list = 0; iline_list < nline_list; iline_list ++){    
        for(long iline = 0; iline < nline_arr[iline_list]; iline ++){
            hd2d_arr[iline_list]->Fill(log10(mu_arr[iline_list][iline]),
                                       beta_arr[iline_list][iline],
                                       ave_arr[iline_list][iline]);
            hd2d_stddev_arr[iline_list]->Fill(log10(mu_arr[iline_list][iline]),
                                              beta_arr[iline_list][iline],
                                              stddev_arr[iline_list][iline]);
            hd2d_index_arr[iline_list]->Fill(imu_arr[iline_list][iline],
                                             ibeta_arr[iline_list][iline],
                                             ave_arr[iline_list][iline]);
            hd2d_stddev_index_arr[iline_list]->Fill(imu_arr[iline_list][iline],
                                                    ibeta_arr[iline_list][iline],
                                                    stddev_arr[iline_list][iline]);
        }
    }

    HistDataNerr2d* hd2d_amean = new HistDataNerr2d;
    HistDataNerr2d* hd2d_stddev = new HistDataNerr2d;
    HistData2dOpe::GetAMean(hd2d_arr, nline_list, hd2d_amean);
    // HistData2dOpe::GetStddev(hd2d_arr, nline_list, hd2d_stddev);
    HistData2dOpe::GetSqrtOfUnbiasedVariance(hd2d_arr, nline_list, hd2d_stddev);

    long ibin_min = hd2d_amean->GetOvalArr()->GetLocValMin();
    printf("(x, y)_min = (%e, %e)\n",
           hd2d_amean->GetHi2d()->GetBinCenterXFromIbin(ibin_min),
           hd2d_amean->GetHi2d()->GetBinCenterYFromIbin(ibin_min));

    // amean of stddev
    double amean_of_stddev = MirMath::GetAMean(hd2d_stddev->GetOvalArr()->GetNdata(),
                                               hd2d_stddev->GetOvalArr()->GetVal());
    printf("amean of stddev = %e\n", amean_of_stddev);


    char outfile[kLineSize];
    sprintf(outfile, "%s/%s.dat",
            argval->GetOutdir().c_str(),
            argval->GetOutfileHead().c_str());
    FILE* fp = fopen(outfile, "w");
    int ndata = hd2d_amean->GetOvalArr()->GetNdata();
    for(int idata = 0; idata < ndata; idata ++){
        fprintf(fp, "%d %e %e\n", idata,
                hd2d_amean->GetOvalArr()->GetValElm(idata),
                hd2d_stddev->GetOvalArr()->GetValElm(idata));
    }
    fclose(fp);
    
    
    double zrange_lo = argval->GetZrangeLo();
    double zrange_up = argval->GetZrangeUp();
    
    TH2D* th2d = hd2d_amean->GenTH2D(0, 0, 0);
    th2d->SetAxisRange(zrange_lo, zrange_up, "Z");
    gStyle->SetPalette(53);
    // TColor::InvertPalette();
    th2d->Draw("COLZ");
    
    
    gPad->Update();
    // TPaletteAxis* palette = (TPaletteAxis*) th2d->GetListOfFunctions()->FindObject("palette");
    // palette->SetX1NDC(0.86);
    // palette->SetX2NDC(0.89);
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

//    /// stddev
//    TH2D* th2d_stddev = hd2d_stddev->GenTH2D(0, 0, 0);
//    th2d_stddev->Draw("COLZ");
//    gPad->Update();
//    th2d_stddev->GetXaxis()->SetTitleSize(0.05);
//    th2d_stddev->GetYaxis()->SetTitleSize(0.05);
//    th2d_stddev->GetXaxis()->SetLabelSize(0.05);
//    th2d_stddev->GetYaxis()->SetLabelSize(0.05);
//
//    char outfig_stddev[kLineSize];
//    sprintf(outfig_stddev, "%s/%s_stddev.png",
//            argval->GetOutdir().c_str(),
//            argval->GetOutfileHead().c_str());
//    
//    root_tool->GetTCanvas()->Print(outfig_stddev);

    delete argval;
    
    return status;
}


