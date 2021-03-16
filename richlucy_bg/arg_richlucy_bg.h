#ifndef MORIIISM_SRT_RICHLUCY_BG_ARG_RICHLUCY_BG_H_
#define MORIIISM_SRT_RICHLUCY_BG_ARG_RICHLUCY_BG_H_

#include "mi_base.h"

class ArgValRichlucyBg : public MiArgBase{
public:
    ArgValRichlucyBg() :
        MiArgBase(),
        progname_(""),
        respdir_(""),
        datafile_(""),
        skyfile_(""),
	bgfile_(""),
        outdir_(""),
        outfile_head_(""),
        nloop_main_(0),
        nloop_em_(0),
        nloop_newton_(0),
        tol_main_(0.0),
        tol_em_(0.0),
        tol_newton_(0.0),
        epsilon_(0.0) {}
    ~ArgValRichlucyBg(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetRespdir() const {return respdir_;};
    string GetDatafile() const {return datafile_;};
    string GetSkyfile() const {return skyfile_;};
    string GetBgfile() const {return bgfile_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};
    int    GetNloopMain() const {return nloop_main_;};
    int    GetNloopEm() const {return nloop_em_;};
    int    GetNloopNewton() const {return nloop_newton_;};
    double GetTolMain() const {return tol_main_;};
    double GetTolEm() const {return tol_em_;};
    double GetTolNewton() const {return tol_newton_;};
    double GetEpsilon() const {return epsilon_;};

private:
    string progname_;
    string respdir_;
    string datafile_;
    string skyfile_;
    string bgfile_;
    string outdir_;
    string outfile_head_;
    int    nloop_main_;
    int    nloop_em_;
    int    nloop_newton_;
    double tol_main_;
    double tol_em_;
    double tol_newton_;
    double epsilon_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_RICHLUCY_BG_ARG_RICHLUCY_BG_H_
