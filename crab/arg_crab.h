#ifndef MORIIISM_SRT_CRAB_ARG_CRAB_H_
#define MORIIISM_SRT_CRAB_ARG_CRAB_H_

#include "mi_base.h"

class ArgValCrab : public MiArgBase{
public:
    ArgValCrab() :
        MiArgBase(),
        progname_(""),
        datafile_(""),
        crabfile_(""),
        skyfile_(""),        
        resp_file_(""),
        eff_file_(""),
        nskyx_(0),
        nskyy_(0),
        ndetx_(0),
        ndety_(0),
        outdir_(""),
        outfile_head_(""),
        nem_(0),
        tol_em_(0.0),
        npm_(0),
        tol_pm_(0.0),
        nnewton_(0),
        tol_newton_(0.0),
        mu_(0.0) {}
    ~ArgValCrab(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetDatafile() const {return datafile_;};
    string GetCrabfile() const {return crabfile_;};
    string GetSkyfile() const {return skyfile_;};    
    string GetRespFile() const {return resp_file_;};
    string GetEffFile() const {return eff_file_;};
    int    GetNskyx() const {return nskyx_;};
    int    GetNskyy() const {return nskyy_;};
    int    GetNdetx() const {return ndetx_;};
    int    GetNdety() const {return ndety_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};
    int    GetNem() const {return nem_;};
    double GetTolEm() const {return tol_em_;};
    int    GetNpm() const {return npm_;};
    double GetTolPm() const {return tol_pm_;};
    int    GetNnewton() const {return nnewton_;};
    double GetTolNewton() const {return tol_newton_;};
    double GetMu() const {return mu_;};
    
private:
    string progname_;
    string datafile_;
    string crabfile_;
    string skyfile_;
    string resp_file_;
    string eff_file_;
    int nskyx_;
    int nskyy_;
    int ndetx_;
    int ndety_;
    string outdir_;
    string outfile_head_;
    int    nem_;
    double tol_em_;
    int    npm_;
    double tol_pm_;    
    int    nnewton_;
    double tol_newton_;
    double mu_;
    
    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_CRAB_ARG_CRAB_H_

