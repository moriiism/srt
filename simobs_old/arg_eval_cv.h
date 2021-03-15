#ifndef MORIIISM_SRT_EVAL_ARG_EVAL_CV_H_
#define MORIIISM_SRT_EVAL_ARG_EVAL_CV_H_

#include "mi_base.h"

class ArgValEvalCv : public MiArgBase{
public:
    ArgValEvalCv() :
        MiArgBase(),
        progname_(""),
        respdir_(""),
        cvlist_(""),
        outdir_(""),
        outfile_head_("") {}
    ~ArgValEvalCv(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetRespdir() const {return respdir_;};
    string GetCvlist() const {return cvlist_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};

private:
    string progname_;
    string respdir_;
    string cvlist_;
    string outdir_;
    string outfile_head_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_EVAL_ARG_EVAL_CV_H_
