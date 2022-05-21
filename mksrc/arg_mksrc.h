#ifndef MORIIISM_SRT_MKSRC_ARG_MKSRC_H_
#define MORIIISM_SRT_MKSRC_ARG_MKSRC_H_

#include "mi_base.h"

class ArgValMksrc : public MiArgBase{
public:
    ArgValMksrc() :
        MiArgBase(),
        progname_(""),
        infile_(""),
        outfile_(""),
        nskyx_(0),
        nskyy_(0){}
    ~ArgValMksrc(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetInfile() const {return infile_;};
    string GetOutfile() const {return outfile_;};
    int    GetNskyx() const {return nskyx_;};
    int    GetNskyy() const {return nskyy_;};    

private:
    string progname_;
    string infile_;
    string outfile_;
    int    nskyx_;
    int    nskyy_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_MKSRC_ARG_MKSRC_H_
