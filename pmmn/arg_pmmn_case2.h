#ifndef MORIIISM_SRT_PMMN_ARG_PMMN_CASE2_H_
#define MORIIISM_SRT_PMMN_ARG_PMMN_CASE2_H_

#include "mi_base.h"

class ArgValPmmnCase2 : public MiArgBase{
public:
    ArgValPmmnCase2() :
        MiArgBase(),
        progname_(""),
        respdir_(""),
        datafile_(""),
        mu_(0.0),
        beta_(0.0),
        outdir_(""),
        outfile_head_(""),
        tol_(0.0),
        nstep_(0),
        lconst_(0.0),
        epsilon_(0.0) {}
    ~ArgValPmmnCase2(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetRespdir() const {return respdir_;};
    string GetDatafile() const {return datafile_;};
    double GetMu() const {return mu_;};
    double GetBeta() const {return beta_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};
    double GetTol() const {return tol_;};
    int    GetNstep() const {return nstep_;};
    double GetLconst() const {return lconst_;};
    double GetEpsilon() const {return epsilon_;};

private:
    string progname_;
    string respdir_;
    string datafile_;
    double mu_;
    double beta_;
    string outdir_;
    string outfile_head_;
    double tol_;
    int    nstep_;
    double lconst_;
    double epsilon_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_PMMN_ARG_PMMN_CASE2_H_
