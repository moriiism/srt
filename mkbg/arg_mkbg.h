#ifndef MORIIISM_SRT_MKBG_ARG_MKBG_H_
#define MORIIISM_SRT_MKBG_ARG_MKBG_H_

#include "mi_base.h"

class ArgValMkbg : public MiArgBase{
public:
    ArgValMkbg() :
        MiArgBase(),
        progname_(""),
        model_file_(""),
        outfile_("") {}
    ~ArgValMkbg(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetModelFile() const {return model_file_;};
    string GetOutfile() const {return outfile_;};

private:
    string progname_;
    string model_file_;
    string outfile_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_MKBG_ARG_MKBG_H_
