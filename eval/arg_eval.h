#ifndef MORIIISM_SRT_EVAL_ARG_EVAL_H_
#define MORIIISM_SRT_EVAL_ARG_EVAL_H_

#include "mi_base.h"

class ArgValEval : public MiArgBase{
public:
    ArgValEval() :
        MiArgBase(),
        progname_(""),
        respdir_(""),
        recfile_(""),
        valfile_(""),
        outdir_(""),
        outfile_head_("") {}
    ~ArgValEval(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetRespdir() const {return respdir_;};
    string GetRecfile() const {return recfile_;};
    string GetValfile() const {return valfile_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};

private:
    string progname_;
    string respdir_;
    string recfile_;
    string valfile_;
    string outdir_;
    string outfile_head_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_EVAL_ARG_EVAL_H_
