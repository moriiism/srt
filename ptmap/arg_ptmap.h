#ifndef MORIIISM_SRT_PTMAP_ARG_PTMAP_H_
#define MORIIISM_SRT_PTMAP_ARG_PTMAP_H_

#include "mi_base.h"

class ArgValPtmap : public MiArgBase{
public:
    ArgValPtmap() :
        MiArgBase(),
        progname_(""),
        cvsmr_list_(""),
        hist_info_file_(""),
        hist_info_index_file_(""),
        zrange_lo_(0.0),
        zrange_up_(0.0),
        outdir_(""),
        outfile_head_("") {}
    ~ArgValPtmap(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetCvsmrList()   const {return cvsmr_list_;};
    string GetHistInfoFile()   const {return hist_info_file_;};
    string GetHistInfoIndexFile()   const {return hist_info_index_file_;};
    double GetZrangeLo() const {return zrange_lo_;};
    double GetZrangeUp() const {return zrange_up_;};
    string GetOutdir() const {return outdir_;};
    string GetOutfileHead()   const {return outfile_head_;};

private:
    string progname_;
    string cvsmr_list_;
    string hist_info_file_;
    string hist_info_index_file_;
    double zrange_lo_;
    double zrange_up_;
    string outdir_;
    string outfile_head_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);    
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_PTMAP_ARG_PTMAP_H_
