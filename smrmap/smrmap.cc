#include "mi_iolib.h"
#include "mir_root_tool.h"
#include "mir_hist2d_nerr.h"
#include "mir_qdp_tool.h"
#include "arg_smrmap.h"
#include "TColor.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[]){
    int status = kRetNormal;
  
    ArgValSmrmap* argval = new ArgValSmrmap;
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
    int* imu_arr = new int [nline];
    int* ibeta_arr = new int [nline];
    double* mu_arr = new double [nline];
    double* beta_arr = new double [nline];
    int* iem_arr = new int [nline];
    int* nzero_arr = new int [nline];
    double* logl_arr = new double [nline];
    double* time_arr = new double [nline];
    double* diff_l_var_arr = new double [nline];

   
    for(long iline = 0; iline < nline; iline ++){
        int nsplit = 0;
        string* split_arr = NULL;
        MiStr::GenSplit(lines_arr[iline], &nsplit, &split_arr);
        printf("nsplit = %d\n", nsplit);
        imu_arr[iline]   = atoi(split_arr[0].c_str());
        ibeta_arr[iline] = atoi(split_arr[1].c_str());
        mu_arr[iline]    = atof(split_arr[2].c_str());
        beta_arr[iline]  = atof(split_arr[3].c_str());
        iem_arr[iline]   = atoi(split_arr[4].c_str());
        nzero_arr[iline] = atoi(split_arr[5].c_str());
        logl_arr[iline]  = atof(split_arr[8].c_str());
        time_arr[iline]  = atof(split_arr[11].c_str());
        diff_l_var_arr[iline] = atof(split_arr[13].c_str());
        MiStr::DelSplit(split_arr);
    }
    
    HistInfo2d* hi2d = new HistInfo2d;
    hi2d->Load(argval->GetHistInfoFile());
    HistDataNerr2d* hd2d_iem = new HistDataNerr2d;
    hd2d_iem->Init(hi2d);
    HistDataNerr2d* hd2d_nzero = new HistDataNerr2d;
    hd2d_nzero->Init(hi2d);
    HistDataNerr2d* hd2d_logl = new HistDataNerr2d;
    hd2d_logl->Init(hi2d);
    HistDataNerr2d* hd2d_time = new HistDataNerr2d;
    hd2d_time->Init(hi2d);
    HistDataNerr2d* hd2d_diff_l_var = new HistDataNerr2d;
    hd2d_diff_l_var->Init(hi2d);
    
    for(long iline = 0; iline < nline; iline ++){
        hd2d_iem->Fill(imu_arr[iline], ibeta_arr[iline], iem_arr[iline]);
        hd2d_nzero->Fill(imu_arr[iline], ibeta_arr[iline], nzero_arr[iline]);
        hd2d_logl->Fill(imu_arr[iline], ibeta_arr[iline], logl_arr[iline]);
        hd2d_time->Fill(imu_arr[iline], ibeta_arr[iline], time_arr[iline]);
        hd2d_diff_l_var->Fill(imu_arr[iline], ibeta_arr[iline], diff_l_var_arr[iline]);
    }

    {
        TH2D* th2d = hd2d_iem->GenTH2D(0, 0, 0);
        // th2d->SetAxisRange(zrange_lo, zrange_up, "Z");
        gStyle->SetPalette(53);
        // TColor::InvertPalette();
        th2d->Draw("COLZ");
    
        gPad->Update();
        th2d->GetXaxis()->SetTitleSize(0.05);
        th2d->GetYaxis()->SetTitleSize(0.05);
        th2d->GetXaxis()->SetLabelSize(0.05);
        th2d->GetYaxis()->SetLabelSize(0.05);

        char outfig[kLineSize];
        sprintf(outfig, "%s/%s_iem.png",
                argval->GetOutdir().c_str(),
                argval->GetOutfileHead().c_str());
        root_tool->GetTCanvas()->Print(outfig);
        delete th2d;
    }

    {
        TH2D* th2d = hd2d_nzero->GenTH2D(0, 0, 0);
        // th2d->SetAxisRange(zrange_lo, zrange_up, "Z");
        gStyle->SetPalette(53);
        // TColor::InvertPalette();
        th2d->Draw("COLZ");
    
        gPad->Update();
        th2d->GetXaxis()->SetTitleSize(0.05);
        th2d->GetYaxis()->SetTitleSize(0.05);
        th2d->GetXaxis()->SetLabelSize(0.05);
        th2d->GetYaxis()->SetLabelSize(0.05);

        char outfig[kLineSize];
        sprintf(outfig, "%s/%s_nzero.png",
                argval->GetOutdir().c_str(),
                argval->GetOutfileHead().c_str());
        root_tool->GetTCanvas()->Print(outfig);
        delete th2d;
    }

    {
        TH2D* th2d = hd2d_logl->GenTH2D(0, 0, 0);
        // th2d->SetAxisRange(zrange_lo, zrange_up, "Z");
        gStyle->SetPalette(53);
        // TColor::InvertPalette();
        th2d->Draw("COLZ");
    
        gPad->Update();
        th2d->GetXaxis()->SetTitleSize(0.05);
        th2d->GetYaxis()->SetTitleSize(0.05);
        th2d->GetXaxis()->SetLabelSize(0.05);
        th2d->GetYaxis()->SetLabelSize(0.05);

        char outfig[kLineSize];
        sprintf(outfig, "%s/%s_logl.png",
                argval->GetOutdir().c_str(),
                argval->GetOutfileHead().c_str());
        root_tool->GetTCanvas()->Print(outfig);
        delete th2d;
    }
    
    {
        TH2D* th2d = hd2d_time->GenTH2D(0, 0, 0);
        // th2d->SetAxisRange(zrange_lo, zrange_up, "Z");
        gStyle->SetPalette(53);
        // TColor::InvertPalette();
        th2d->Draw("COLZ");
    
        gPad->Update();
        th2d->GetXaxis()->SetTitleSize(0.05);
        th2d->GetYaxis()->SetTitleSize(0.05);
        th2d->GetXaxis()->SetLabelSize(0.05);
        th2d->GetYaxis()->SetLabelSize(0.05);

        char outfig[kLineSize];
        sprintf(outfig, "%s/%s_time.png",
                argval->GetOutdir().c_str(),
                argval->GetOutfileHead().c_str());
        root_tool->GetTCanvas()->Print(outfig);
        delete th2d;
    }

    {
        TH2D* th2d = hd2d_diff_l_var->GenTH2D(0, 0, 0);
        // th2d->SetAxisRange(zrange_lo, zrange_up, "Z");
        gStyle->SetPalette(53);
        // TColor::InvertPalette();
        th2d->Draw("COLZ");
    
        gPad->Update();
        th2d->GetXaxis()->SetTitleSize(0.05);
        th2d->GetYaxis()->SetTitleSize(0.05);
        th2d->GetXaxis()->SetLabelSize(0.05);
        th2d->GetYaxis()->SetLabelSize(0.05);

        char outfig[kLineSize];
        sprintf(outfig, "%s/%s_diff_l_var.png",
                argval->GetOutdir().c_str(),
                argval->GetOutfileHead().c_str());
        root_tool->GetTCanvas()->Print(outfig);
        delete th2d;
    }
    
    
    
    delete argval;
    
    return status;
}


