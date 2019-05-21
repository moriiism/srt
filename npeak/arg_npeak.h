#ifndef MORIIISM_SRT_NPEAK_ARG_NPEAK_H_
#define MORIIISM_SRT_NPEAK_ARG_NPEAK_H_

#include "mi_base.h"

// xpos, ypos: origin of the new x-y frame.
// theta     : angle of the new x axis direction
//             from original x axis (deg)
// y_lo, y_up: range of the new y axis 

class ArgValNpeak : public MiArgBase{
public:
    ArgValNpeak() :
        MiArgBase(),
        progname_(""),
        infile_(""),
        xpos_(0.0),
        ypos_(0.0),
        theta_(0.0),
        y_lo_(0.0),
        y_up_(0.0),
        nbinx_new_(0),
        xlo_new_(0.0),
        xup_new_(0.0),
        significance_(0.0),
        outdir_(""),
        outfile_head_("") {}
    ~ArgValNpeak(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetInfile() const {return infile_;};
    double GetXpos() const {return xpos_;};
    double GetYpos() const {return ypos_;};
    double GetTheta() const {return theta_;};
    double GetYLo() const {return y_lo_;};
    double GetYUp() const {return y_up_;};
    int    GetNbinxNew() const {return nbinx_new_;};
    double GetXloNew() const {return xlo_new_;};
    double GetXupNew() const {return xup_new_;};
    double GetSignificance() const {return significance_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead() const {return outfile_head_;};

private:
    string progname_;
    string infile_;
    double xpos_;
    double ypos_;
    double theta_;
    double y_lo_;
    double y_up_;
    int    nbinx_new_;
    double xlo_new_;
    double xup_new_;
    double significance_;
    string outdir_;
    string outfile_head_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_NPEAK_ARG_NPEAK_H_
