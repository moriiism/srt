#ifndef MORIIISM_SRT_SIMOBS_ARG_SIMOBS_H_
#define MORIIISM_SRT_SIMOBS_ARG_SIMOBS_H_

#include "mi_base.h"

class ArgValSimobs : public MiArgBase{
public:
    ArgValSimobs() :
        MiArgBase(),
        progname_(""),
        respdir_(""),
        infile_(""),
        nevt_(0),
        rand_seed_sky_(0),
        rand_seed_det_(0),
        rand_seed_partition_(0),
        nfold_(0),
        npartition_(0),
        outdir_(""),
        outfile_head_("") {}
    ~ArgValSimobs(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetRespdir() const {return respdir_;};
    string GetInfile() const {return infile_;};
    int    GetNevt() const {return nevt_;};
    int    GetRandSeedSky() const {return rand_seed_sky_;};
    int    GetRandSeedDet() const {return rand_seed_det_;};
    int    GetRandSeedPartition() const {return rand_seed_partition_;};
    int    GetNfold() const {return nfold_;};
    int    GetNpartition() const {return npartition_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};

private:
    string progname_;
    string respdir_;
    string infile_;
    int    nevt_;
    int    rand_seed_sky_;
    int    rand_seed_det_;
    int    rand_seed_partition_;
    int    nfold_;
    int    npartition_;
    string outdir_;
    string outfile_head_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_SIMOBS_ARG_SIMOBS_H_
