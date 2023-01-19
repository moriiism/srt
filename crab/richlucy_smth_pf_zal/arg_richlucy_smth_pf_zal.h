#ifndef MORIIISM_SRT_CRAB_RICHLUCY_SMTH_PF_ZAL_ARG_RICHLUCY_SMTH_PF_ZAL_H_
#define MORIIISM_SRT_CRAB_RICHLUCY_SMTH_PF_ZAL_ARG_RICHLUCY_SMTH_PF_ZAL_H_

#include "mi_base.h"

class ArgValRichlucySmthPfZal : public MiArgBase{
public:
    ArgValRichlucySmthPfZal() :
        MiArgBase(),
        progname_(""),
        data_list_(""),
        bg_file_(""),
        fixed_src_norm_file_(""),
        resp_norm_file_(""),
        eff_file_(""),
        nskyx_(0),
        nskyy_(0),
        ndetx_(0),
        ndety_(0),
        outdir_(""),
        outfile_head_(""),
        nem_(0),
        tol_em_(0.0),
        mu_(0.0),
        gamma_(0.0),
        acc_method_(""),
        nfold_cv_(0){}
    ~ArgValRichlucySmthPfZal(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetDataList() const {return data_list_;};
    string GetBgFile() const {return bg_file_;};
    string GetFixedSrcNormFile() const {return fixed_src_norm_file_;};
    string GetRespNormFile() const {return resp_norm_file_;};
    string GetEffFile() const {return eff_file_;};
    int    GetNskyx() const {return nskyx_;};
    int    GetNskyy() const {return nskyy_;};
    int    GetNdetx() const {return ndetx_;};
    int    GetNdety() const {return ndety_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};
    int    GetNem() const {return nem_;};
    double GetTolEm() const {return tol_em_;};
    double GetMu() const {return mu_;};
    double GetGamma() const {return gamma_;};
    string GetAccMethod() const {return acc_method_;};
    int    GetNfoldCv() const {return nfold_cv_;};

private:
    string progname_;
    string data_list_;
    string bg_file_;
    string fixed_src_norm_file_;
    string resp_norm_file_;
    string eff_file_;
    int nskyx_;
    int nskyy_;
    int ndetx_;
    int ndety_;
    string outdir_;
    string outfile_head_;
    int    nem_;
    double tol_em_;
    double mu_;
    double gamma_;
    string acc_method_;
    int nfold_cv_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_CRAB_RICHLUCY_SMTH_PF_ZAL_ARG_RICHLUCY_SMTH_PF_ZAL_H_
