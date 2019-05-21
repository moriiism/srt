#ifndef MORIIISM_SRT_RICHLUCY_ARG_RICHLUCY_H_
#define MORIIISM_SRT_RICHLUCY_ARG_RICHLUCY_H_

#include "mi_base.h"

class ArgValRichlucy : public MiArgBase{
public:
    ArgValRichlucy() :
        MiArgBase(),
        progname_(""),
        respdir_(""),
        datafile_(""),
        skyfile_(""),
        outdir_(""),
        outfile_head_(""),
        nem_(0),
        tol_em_(0.0),
        tol_diff_l_var_(0.0),
        flag_line_search_(0),
        epsilon_(0.0) {}
    ~ArgValRichlucy(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetRespdir() const {return respdir_;};
    string GetDatafile() const {return datafile_;};
    string GetSkyfile() const {return skyfile_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};
    int    GetNem() const {return nem_;};
    double GetTolEm() const {return tol_em_;};
    double GetTolDiffLVar() const {return tol_diff_l_var_;};
    int    GetFlagLineSearch() const {return flag_line_search_;};
    double GetEpsilon() const {return epsilon_;};

private:
    string progname_;
    string respdir_;
    string datafile_;
    string skyfile_;
    string outdir_;
    string outfile_head_;
    int    nem_;
    double tol_em_;
    double tol_diff_l_var_;
    int    flag_line_search_;
    double epsilon_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_RICHLUCY_ARG_RICHLUCY_H_
