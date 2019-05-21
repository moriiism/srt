#ifndef MORIIISM_SRT_HPD_ARG_HPD_H_
#define MORIIISM_SRT_HPD_ARG_HPD_H_

#include "mi_base.h"

class ArgValHpd : public MiArgBase{
public:
    ArgValHpd() :
        MiArgBase(),
        progname_(""),
        infile_(""),
        xpos_(0.0),
        ypos_(0.0),
        outfile_("") {}
    ~ArgValHpd(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetInfile() const {return infile_;};
    double GetXpos() const {return xpos_;};
    double GetYpos() const {return ypos_;};
    string GetOutfile() const {return outfile_;};

private:
    string progname_;
    string infile_;
    double xpos_;
    double ypos_;
    string outfile_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_HPD_ARG_HPD_H_
