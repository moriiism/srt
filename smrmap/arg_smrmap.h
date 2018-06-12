#ifndef MORIIISM_SRT_SMRMAP_ARG_SMRMAP_H_
#define MORIIISM_SRT_SMRMAP_ARG_SMRMAP_H_

#include "mi_base.h"

class ArgValSmrmap : public MiArgBase{
public:
    ArgValSmrmap() :
        MiArgBase(),
        progname_(""),
        infile_(""),
        hist_info_file_(""),
        outdir_(""),
        outfile_head_("") {}
    ~ArgValSmrmap(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetInfile()   const {return infile_;};
    string GetHistInfoFile()   const {return hist_info_file_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead()   const {return outfile_head_;};

private:
    string progname_;
    string infile_;
    string hist_info_file_;
    string outdir_;
    string outfile_head_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);    
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_SMRMAP_ARG_SMRMAP_H_
