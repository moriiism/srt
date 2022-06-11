#ifndef MORIIISM_SRT_MKOBS_CV_ARG_MKOBS_CV_H_
#define MORIIISM_SRT_MKOBS_CV_ARG_MKOBS_CV_H_

#include "mi_base.h"

class ArgValMkobsCv : public MiArgBase{
public:
    ArgValMkobsCv() :
        MiArgBase(),
        progname_(""),
        orgfile_(""),
        rand_seed_(0),
        nfold_(0),
        outdir_(""),
        outfile_head_("") {}
    ~ArgValMkobsCv(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetOrgfile() const {return orgfile_;};
    int    GetRandSeed() const {return rand_seed_;};
    int    GetNfold() const {return nfold_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};

private:
    string progname_;
    string orgfile_;
    int    rand_seed_;
    int    nfold_;
    string outdir_;
    string outfile_head_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_MKOBS_CV_ARG_MKOBS_CV_H_
