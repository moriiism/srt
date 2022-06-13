#ifndef MORIIISM_SRT_RICHLUCY_CRAB_ARG_EVAL_VAL_H_
#define MORIIISM_SRT_RICHLUCY_CRAB_ARG_EVAL_VAL_H_

#include "mi_base.h"

class ArgValEvalVal : public MiArgBase{
public:
    ArgValEvalVal() :
        MiArgBase(),
        progname_(""),
        resp_file_(""),
        recfile_(""),
        valfile_(""),
        outdir_(""),
        outfile_head_("") {}
    ~ArgValEvalVal(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetRespFile() const {return resp_file_;};
    string GetRecfile() const {return recfile_;};
    string GetValfile() const {return valfile_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};

private:
    string progname_;
    string resp_file_;
    string recfile_;
    string valfile_;
    string outdir_;
    string outfile_head_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_RICHLUCY_CRAB_ARG_EVAL_VAL_H_
