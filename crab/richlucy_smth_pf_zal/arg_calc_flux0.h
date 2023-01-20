#ifndef MORIIISM_SRT_CRAB_RICHLUCY_SMTH_PF_ZAL_ARG_CALC_FLUX0_H_
#define MORIIISM_SRT_CRAB_RICHLUCY_SMTH_PF_ZAL_ARG_CALC_FLUX0_H_

#include "mi_base.h"

class ArgValCalcFlux0 : public MiArgBase{
public:
    ArgValCalcFlux0() :
        MiArgBase(),
        progname_(""),
        data_on_file_(""),
        data_off_file_(""),
        phase_ratio_on_(0.0),
        phase_ratio_off_(0.0),
        live_time_ratio_on_(0.0),
        live_time_ratio_off_(0.0) {}
    ~ArgValCalcFlux0(){
        Null();
    }
    void Init(int argc, char* argv[]);
    void Print(FILE* fp) const;

    string GetProgname() const {return progname_;};
    string GetDataOnFile() const {return data_on_file_;};
    string GetDataOffFile() const {return data_off_file_;};
    double GetPhaseRatioOn() const {return phase_ratio_on_;};
    double GetPhaseRatioOff() const {return phase_ratio_off_;};
    double GetLiveTimeRatioOn() const {
        return live_time_ratio_on_;};
    double GetLiveTimeRatioOff() const {
        return live_time_ratio_off_;};

private:
    string progname_;
    string data_on_file_;
    string data_off_file_;
    double phase_ratio_on_;
    double phase_ratio_off_;
    double live_time_ratio_on_;
    double live_time_ratio_off_;

    void Null();
    void SetOption(int argc, char* argv[], option* long_options);
    void Usage(FILE* fp) const;
};

#endif // MORIIISM_SRT_CRAB_RICHLUCY_SMTH_PF_ZAL_ARG_CALC_FLUX0_H_
