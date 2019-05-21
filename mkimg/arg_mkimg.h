#ifndef MORIIISM_SRT_MKIMG_ARG_MKIMG_H_
#define MORIIISM_SRT_MKIMG_ARG_MKIMG_H_

#include "mi_base.h"

class ArgValMkimg : public MiArgBase{
public:
    ArgValMkimg() :
        MiArgBase(),
        progname_(""),
        infile_(""),
        dotfile_(""),
        outdir_(""),
        outfile_head_("") {}
    ~ArgValMkimg(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetInfile() const {return infile_;};
    string GetDotfile() const {return dotfile_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};

private:
    string progname_;
    string infile_;
    string dotfile_;
    string outdir_;
    string outfile_head_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_MKIMG_ARG_MKIMG_H_
