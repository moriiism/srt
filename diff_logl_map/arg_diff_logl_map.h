#ifndef MORIIISM_SRT_DIFF_LOGL_MAP_ARG_DIFF_LOGL_MAP_H_
#define MORIIISM_SRT_DIFF_LOGL_MAP_ARG_DIFF_LOGL_MAP_H_

#include "mi_base.h"

class ArgValDiffLoglMap : public MiArgBase{
public:
    ArgValDiffLoglMap() :
        MiArgBase(),
        progname_(""),
        infile_flat_(""),
        infile_init_(""),
        hist_info_file_(""),
        outdir_(""),
        outfile_head_("") {}
    ~ArgValDiffLoglMap(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetInfileFlat()   const {return infile_flat_;};
    string GetInfileInit()   const {return infile_init_;};
    string GetHistInfoFile()   const {return hist_info_file_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead()   const {return outfile_head_;};

private:
    string progname_;
    string infile_flat_;
    string infile_init_;
    string hist_info_file_;
    string outdir_;
    string outfile_head_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);    
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_DIFF_LOGL_MAP_ARG_DIFF_LOGL_MAP_H_
