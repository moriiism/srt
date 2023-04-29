#ifndef MORIIISM_SRT_MKOBS_NPH_ARG_MKOBS_NPH_RATIO_H_
#define MORIIISM_SRT_MKOBS_NPH_ARG_MKOBS_NPH_RATIO_H_

#include "mi_base.h"

class ArgValMkobsNphRatio: public MiArgBase{
public:
    ArgValMkobsNphRatio() :
        MiArgBase(),
        progname_(""),
        orgfile_(""),
        rand_seed_(0),
        ratio_(0.0),
        outdir_(""),
        outfile_head_("") {}
    ~ArgValMkobsNphRatio(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetOrgfile() const {return orgfile_;};
    int    GetRandSeed() const {return rand_seed_;};
    double GetRatio() const {return ratio_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};

private:
    string progname_;
    string orgfile_;
    int    rand_seed_;
    double ratio_;
    string outdir_;
    string outfile_head_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_MKOBS_NPH_ARG_MKOBS_NPH_RATIO_H_
