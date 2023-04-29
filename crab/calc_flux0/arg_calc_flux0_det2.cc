#include "arg_calc_flux0_det2.h"

// public

void ArgValCalcFlux0Det2::Init(int argc, char* argv[])
{
    progname_ = "calc_flux0_det2";

    option long_options[] = {
        {"debug",       required_argument, NULL, 'd'},
        {"help",        required_argument, NULL, 'h'},
        {"verbose",     required_argument, NULL, 'v'},
        {0, 0, 0, 0}
    };

    // long option default

    SetOption(argc, argv, long_options);
    
    printf("ArgVal::Init: # of arg = %d\n", argc - optind);
    int narg = 10;
    if (argc - optind != narg){
        printf("# of arguments must be %d.\n", narg);
        Usage(stdout);
    }
    int iarg = optind;
    data_on_file_det1_        = argv[iarg]; iarg++;
    data_on_file_det2_        = argv[iarg]; iarg++;
    live_time_ratio_on_det1_  = atof(argv[iarg]); iarg++;
    live_time_ratio_on_det2_  = atof(argv[iarg]); iarg++;
    phase_ratio_on_           = atof(argv[iarg]); iarg++;
    data_off_file_det1_       = argv[iarg]; iarg++;
    data_off_file_det2_       = argv[iarg]; iarg++;
    live_time_ratio_off_det1_ = atof(argv[iarg]); iarg++;
    live_time_ratio_off_det2_ = atof(argv[iarg]); iarg++;
    phase_ratio_off_          = atof(argv[iarg]); iarg++;
}

void ArgValCalcFlux0Det2::Print(FILE* fp) const
{
    fprintf(fp, "%s: g_flag_debug   : %d\n", __func__,
            g_flag_debug);
    fprintf(fp, "%s: g_flag_help    : %d\n", __func__,
            g_flag_help);
    fprintf(fp, "%s: g_flag_verbose : %d\n", __func__,
            g_flag_verbose);

    fprintf(fp, "%s: progname_       : %s\n", __func__,
            progname_.c_str());
    fprintf(fp, "%s: data_on_file_det1_ : %s\n", __func__,
            data_on_file_det1_.c_str());
    fprintf(fp, "%s: data_on_file_det2_ : %s\n", __func__,
            data_on_file_det2_.c_str());
    fprintf(fp, "%s: live_time_ratio_on_det1_ : %e\n", __func__,
            live_time_ratio_on_det1_);
    fprintf(fp, "%s: live_time_ratio_on_det2_ : %e\n", __func__,
            live_time_ratio_on_det2_);
    fprintf(fp, "%s: phase_ratio_on_ : %e\n", __func__,
            phase_ratio_on_);
    fprintf(fp, "%s: data_off_file_det1_ : %s\n", __func__,
            data_off_file_det1_.c_str());
    fprintf(fp, "%s: data_off_file_det2_ : %s\n", __func__,
            data_off_file_det2_.c_str());
    fprintf(fp, "%s: live_time_ratio_off_det1_ : %e\n", __func__,
            live_time_ratio_off_det1_);
    fprintf(fp, "%s: live_time_ratio_off_det2_ : %e\n", __func__,
            live_time_ratio_off_det2_);
    fprintf(fp, "%s: phase_ratio_off_ : %e\n", __func__,
            phase_ratio_off_);
}


// private

void ArgValCalcFlux0Det2::Null()
{
    progname_ = "";
    data_on_file_det1_ = "";
    data_on_file_det2_ = "";
    live_time_ratio_on_det1_ = 0.0;
    live_time_ratio_on_det2_ = 0.0;
    phase_ratio_on_ = 0.0;
    data_off_file_det1_ = "";
    data_off_file_det2_ = "";
    live_time_ratio_off_det1_ = 0.0;
    live_time_ratio_off_det2_ = 0.0;    
    phase_ratio_off_ = 0.0;
}

void ArgValCalcFlux0Det2::SetOption(
    int argc, char* argv[], option* long_options)
{
    if(0 < g_flag_verbose){
        MPrintInfo("start...");
    }
    // option default
    g_flag_debug   = 0;
    g_flag_help    = 0;
    g_flag_verbose = 0;
    while (1) {
        int option_index = 0;
        int retopt = getopt_long(argc, argv, "dhv",
                                 long_options, &option_index);
        if(-1 == retopt)
            break;
        switch (retopt) {
        case 0:
            // long option
            break;
        case 'd':
            g_flag_debug = atoi(optarg);
            printf("%s: g_flag_debug = %d\n", __func__, g_flag_debug);
            break;
        case 'h':
            g_flag_help = atoi(optarg);
            printf("%s: g_flag_help = %d\n", __func__, g_flag_help);
            if(0 != g_flag_help){
                Usage(stdout);
            }
            break;
        case 'v':
            g_flag_verbose = atoi(optarg);
            printf("%s: g_flag_verbose = %d\n", __func__, g_flag_verbose);
            break;
        case '?':
            printf("%s: retopt (= %c) is invalid flag.\n",
                   __func__, retopt);
            Usage(stdout);
            break;
        default:
            printf("%s: error: getopt returned character code 0%o ??\n", __func__, retopt);
            abort();
        }
    }
    if(0 < g_flag_verbose){
        MPrintInfo("done.");
    }
}


void ArgValCalcFlux0Det2::Usage(FILE* fp) const
{
    fprintf(fp,
            "usage: %s [--help (0)] [--verbose (0)] [--debug (0)] "
            "data_on_file_det1  data_on_file_det2  "
            "live_time_ratio_on_det1  live_time_ratio_on_det2  "
            "phase_ratio_on  "
            "data_off_file_det1  data_off_file_det2  "
            "live_time_ratio_off_det1  live_time_ratio_off_det2  "
            "phase_ratio_off\n",
            progname_.c_str());
    abort();
}
