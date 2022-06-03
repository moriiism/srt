#ifndef MORIIISM_SRT_MKIMG_POINTS_ARG_MKIMG_POINTS_H_
#define MORIIISM_SRT_MKIMG_POINTS_ARG_MKIMG_POINTS_H_

#include "mi_base.h"

class ArgValMkimgPoints: public MiArgBase{
public:
    ArgValMkimgPoints():
        MiArgBase(),
        progname_(""),
        infile_(""),
        outdir_(""),
        outfile_head_(""),
        nskyx_(0),
        nskyy_(0){}
    ~ArgValMkimgPoints(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetInfile() const {return infile_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};
    int    GetNskyx() const {return nskyx_;};
    int    GetNskyy() const {return nskyy_;};

private:
    string progname_;
    string infile_;
    string outdir_;
    string outfile_head_;
    int    nskyx_;
    int    nskyy_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_MKIMG_POINTS_ARG_MKIMG_POINTS_H_
