#ifndef MORIIISM_SRT_MKOBS_NPH_ARG_MKOBS_NPH_H_
#define MORIIISM_SRT_MKOBS_NPH_ARG_MKOBS_NPH_H_

#include "mi_base.h"

class ArgValMkobsNph : public MiArgBase{
public:
    ArgValMkobsNph() :
        MiArgBase(),
        progname_(""),
        orgfile_(""),
        rand_seed_(0),
        nphoton_(0),
        outdir_(""),
        outfile_head_("") {}
    ~ArgValMkobsNph(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetOrgfile() const {return orgfile_;};
    int    GetRandSeed() const {return rand_seed_;};
    int    GetNphoton() const {return nphoton_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};

private:
    string progname_;
    string orgfile_;
    int    rand_seed_;
    int    nphoton_;
    string outdir_;
    string outfile_head_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_MKOBS_NPH_ARG_MKOBS_NPH_H_
