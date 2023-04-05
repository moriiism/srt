#ifndef MORIIISM_SRT_MKRESP_DET2_ARG_MKRESP_DET2_H_
#define MORIIISM_SRT_MKRESP_DET2_ARG_MKRESP_DET2_H_

#include "mi_base.h"

class ArgValMkrespDet2 : public MiArgBase{
public:
    ArgValMkrespDet2() :
        MiArgBase(),
        progname_(""),
        respdir1_(""),
        respdir2_(""),
        outdir_(""),
        outfile_head_(""),
        nskyx_(0),
        nskyy_(0),
        nphoton_input_(0){}
    ~ArgValMkrespDet2(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetRespdir1() const {return respdir1_;};
    string GetRespdir2() const {return respdir2_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};
    int    GetNskyx() const {return nskyx_;};
    int    GetNskyy() const {return nskyy_;};
    int    GetNphotonInput() const {return nphoton_input_;};

private:
    string progname_;
    string respdir1_;
    string respdir2_;    
    string outdir_;
    string outfile_head_;
    int    nskyx_;
    int    nskyy_;
    int    nphoton_input_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_MKRESP_DET2_ARG_MKRESP_DET2_H_
