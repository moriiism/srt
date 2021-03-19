#ifndef MORIIISM_SRT_SIMOBS_ARG_SIMOBS_H_
#define MORIIISM_SRT_SIMOBS_ARG_SIMOBS_H_

#include "mi_base.h"

class ArgValSimobs : public MiArgBase{
public:
    ArgValSimobs() :
        MiArgBase(),
        progname_(""),
        respfile_(""),
        srcfile_(""),
        nevt_src_(0),
        bgfile_(""),
        nevt_bg_(0),
        rand_seed_det_(0),
        outdir_(""),
        outfile_head_("") {}
    ~ArgValSimobs(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetRespfile() const {return respfile_;};
    string GetSrcfile() const {return srcfile_;};
    int    GetNevtSrc() const {return nevt_src_;};
    string GetBgfile() const {return bgfile_;};
    int    GetNevtBg() const {return nevt_bg_;};
    int    GetRandSeedDet() const {return rand_seed_det_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};

private:
    string progname_;
    string respfile_;
    string srcfile_;
    int    nevt_src_;
    string bgfile_;
    int    nevt_bg_;
    int    rand_seed_det_;
    string outdir_;
    string outfile_head_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_SIMOBS_ARG_SIMOBS_H_
