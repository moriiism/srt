#ifndef MORIIISM_SRT_RICHLUCY_DET2_ARG_RICHLUCY_DET2_H_
#define MORIIISM_SRT_RICHLUCY_DET2_ARG_RICHLUCY_DET2_H_

#include "mi_base.h"

class ArgValRichlucyDet2 : public MiArgBase{
public:
    ArgValRichlucyDet2() :
        MiArgBase(),
        progname_(""),
        datafile1_(""),
        datafile2_(""),
        bg_file1_(""),
        bg_file2_(""),
        resp_norm_file1_(""),
        resp_norm_file2_(""),
        eff_file_(""),
        nskyx_(0),
        nskyy_(0),
        ndetx_(0),
        ndety_(0),
        outdir_(""),
        outfile_head_(""),
        nem_(0),
        tol_em_(0.0),
        acc_method_("") {}
    ~ArgValRichlucyDet2(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetDatafile1() const {return datafile1_;};
    string GetDatafile2() const {return datafile2_;};
    string GetBgFile1() const {return bg_file1_;};
    string GetBgFile2() const {return bg_file2_;};    
    string GetRespNormFile1() const {return resp_norm_file1_;};
    string GetRespNormFile2() const {return resp_norm_file2_;};    
    string GetEffFile() const {return eff_file_;};
    int    GetNskyx() const {return nskyx_;};
    int    GetNskyy() const {return nskyy_;};
    int    GetNdetx() const {return ndetx_;};
    int    GetNdety() const {return ndety_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};
    int    GetNem() const {return nem_;};
    double GetTolEm() const {return tol_em_;};
    string GetAccMethod() const {return acc_method_;};

private:
    string progname_;
    string datafile1_;
    string datafile2_;
    string bg_file1_;
    string bg_file2_;    
    string resp_norm_file1_;
    string resp_norm_file2_;    
    string eff_file_;
    int nskyx_;
    int nskyy_;
    int ndetx_;
    int ndety_;
    string outdir_;
    string outfile_head_;
    int    nem_;
    double tol_em_;
    string acc_method_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_RICHLUCY_DET2_ARG_RICHLUCY_DET2_H_
