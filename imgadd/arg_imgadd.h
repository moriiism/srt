#ifndef MORIIISM_SRT_IMGADD_ARG_IMGADD_H_
#define MORIIISM_SRT_IMGADD_ARG_IMGADD_H_

#include "mi_base.h"

class ArgValImgadd : public MiArgBase{
public:
    ArgValImgadd() :
        MiArgBase(),
        progname_(""),
        imgfile1_(""),
        imgfile2_(""),
        coeff1_(0.0),
        coeff2_(0.0),
        outfile_("") {}
    ~ArgValImgadd(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetImgfile1() const {return imgfile1_;};
    string GetImgfile2() const {return imgfile2_;};
    double GetCoeff1() const {return coeff1_;};
    double GetCoeff2() const {return coeff2_;};
    string GetOutfile() const {return outfile_;};

private:
    string progname_;
    string imgfile1_;
    string imgfile2_;
    double coeff1_;
    double coeff2_;
    string outfile_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_IMGADD_ARG_IMGADD_H_
