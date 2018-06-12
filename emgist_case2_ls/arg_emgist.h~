#ifndef MORIIISM_SRT_EMGIST_ARG_EMGIST_H_
#define MORIIISM_SRT_EMGIST_ARG_EMGIST_H_

#include "mi_base.h"

class ArgValEmgist : public MiArgBase{
public:
    ArgValEmgist() :
        MiArgBase(),
        progname_(""),
        respdir_(""),
        datafile_(""),
        skyfile_(""),
        mu_(0.0),
        beta_(0.0),
        outdir_(""),
        outfile_head_(""),
        nem_(0),
        tol_em_(0.0),
        npm_(0),
        tol_pm_(0.0),
        flag_line_search_(0),
        lconst_(0.0),
        epsilon_(0.0) {}
    ~ArgValEmgist(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetRespdir() const {return respdir_;};
    string GetDatafile() const {return datafile_;};
    string GetSkyfile() const {return skyfile_;};
    double GetMu() const {return mu_;};
    double GetBeta() const {return beta_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};
    int    GetNem() const {return nem_;};
    double GetTolEm() const {return tol_em_;};
    int    GetNpm() const {return npm_;};
    double GetTolPm() const {return tol_pm_;};
    int    GetFlagLineSearch() const {return flag_line_search_;};
    double GetLconst() const {return lconst_;};
    double GetEpsilon() const {return epsilon_;};

private:
    string progname_;
    string respdir_;
    string datafile_;
    string skyfile_;
    double mu_;
    double beta_;
    string outdir_;
    string outfile_head_;
    int    nem_;
    double tol_em_;
    int    npm_;
    double tol_pm_;    
    int    flag_line_search_;
    double lconst_;
    double epsilon_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_EMGIST_ARG_EMGIST_H_
