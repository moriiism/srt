#ifndef MORIIISM_SRT_MKBG_ARG_MKBG_H_
#define MORIIISM_SRT_MKBG_ARG_MKBG_H_

#include "mi_base.h"

class ArgValMkbg : public MiArgBase{
public:
    ArgValMkbg() :
        MiArgBase(),
        progname_(""),
        model_file_(""),
        outfile_(""),
        ndetx_(0),
        ndety_(0) {}
    ~ArgValMkbg(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetModelFile() const {return model_file_;};
    string GetOutfile() const {return outfile_;};
    int    GetNdetx() const {return ndetx_;};
    int    GetNdety() const {return ndety_;};    

private:
    string progname_;
    string model_file_;
    string outfile_;
    int    ndetx_;
    int    ndety_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_MKBG_ARG_MKBG_H_
