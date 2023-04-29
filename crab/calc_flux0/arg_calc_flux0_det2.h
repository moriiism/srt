#ifndef MORIIISM_SRT_CRAB_CALC_FLUX0_ARG_CALC_FLUX0_DET2_H_
#define MORIIISM_SRT_CRAB_CALC_FLUX0_ARG_CALC_FLUX0_DET2_H_

#include "mi_base.h"

class ArgValCalcFlux0Det2 : public MiArgBase{
public:
    ArgValCalcFlux0Det2() :
        MiArgBase(),
        progname_(""),
        data_on_file_det1_(""),
        data_on_file_det2_(""),
        live_time_ratio_on_det1_(0.0),
        live_time_ratio_on_det2_(0.0),
        phase_ratio_on_(0.0),
        data_off_file_det1_(""),
        data_off_file_det2_(""),
        live_time_ratio_off_det1_(0.0),
        live_time_ratio_off_det2_(0.0),
        phase_ratio_off_(0.0) {}
    ~ArgValCalcFlux0Det2(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetDataOnFileDet1() const {return data_on_file_det1_;};
    string GetDataOnFileDet2() const {return data_on_file_det2_;};
    double GetLiveTimeRatioOnDet1() const {
        return live_time_ratio_on_det1_;};
    double GetLiveTimeRatioOnDet2() const {
        return live_time_ratio_on_det2_;};    
    double GetPhaseRatioOn() const {return phase_ratio_on_;};
    string GetDataOffFileDet1() const {return data_off_file_det1_;};
    string GetDataOffFileDet2() const {return data_off_file_det2_;};
    double GetLiveTimeRatioOffDet1() const {
        return live_time_ratio_off_det1_;};
    double GetLiveTimeRatioOffDet2() const {
        return live_time_ratio_off_det2_;};
    double GetPhaseRatioOff() const {return phase_ratio_off_;};

private:
    string progname_;
    string data_on_file_det1_;
    string data_on_file_det2_;
    double live_time_ratio_on_det1_;
    double live_time_ratio_on_det2_;    
    double phase_ratio_on_;
    string data_off_file_det1_;
    string data_off_file_det2_;
    double live_time_ratio_off_det1_;
    double live_time_ratio_off_det2_;
    double phase_ratio_off_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_CRAB_CALC_FLUX0_ARG_CALC_FLUX0_DET2_H_
