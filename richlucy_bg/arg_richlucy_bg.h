#ifndef MORIIISM_SRT_RICHLUCY_BG_ARG_RICHLUCY_BG_H_
#define MORIIISM_SRT_RICHLUCY_BG_ARG_RICHLUCY_BG_H_

#include "mi_base.h"

class ArgValRichlucyBg : public MiArgBase{
public:
    ArgValRichlucyBg() :
        MiArgBase(),
        progname_(""),
        datafile_(""),
        skyfile_(""),
	bgfile_(""),
        resp_file_(""),
        eff_file_(""),
        nskyx_(0),
        nskyy_(0),
        ndetx_(0),
        ndety_(0),
        outdir_(""),
        outfile_head_(""),
        nloop_main_(0),
        nloop_em_(0),
        nloop_newton_(0),
        tol_main_(0.0),
        tol_em_(0.0),
        tol_newton_(0.0) {}
    ~ArgValRichlucyBg(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetDatafile() const {return datafile_;};
    string GetSkyfile() const {return skyfile_;};
    string GetBgfile() const {return bgfile_;};
    string GetRespFile() const {return resp_file_;};
    string GetEffFile() const {return eff_file_;};
    int    GetNskyx() const {return nskyx_;};
    int    GetNskyy() const {return nskyy_;};
    int    GetNdetx() const {return ndetx_;};
    int    GetNdety() const {return ndety_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};
    int    GetNloopMain() const {return nloop_main_;};
    int    GetNloopEm() const {return nloop_em_;};
    int    GetNloopNewton() const {return nloop_newton_;};
    double GetTolMain() const {return tol_main_;};
    double GetTolEm() const {return tol_em_;};
    double GetTolNewton() const {return tol_newton_;};

private:
    string progname_;
    string datafile_;
    string skyfile_;
    string bgfile_;
    string resp_file_;
    string eff_file_;
    int nskyx_;
    int nskyy_;
    int ndetx_;
    int ndety_;
    string outdir_;
    string outfile_head_;
    int    nloop_main_;
    int    nloop_em_;
    int    nloop_newton_;
    double tol_main_;
    double tol_em_;
    double tol_newton_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_RICHLUCY_BG_ARG_RICHLUCY_BG_H_
