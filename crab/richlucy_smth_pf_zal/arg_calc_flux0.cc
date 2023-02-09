#include "arg_calc_flux0.h"

// public

void ArgValCalcFlux0::Init(int argc, char* argv[])
{
    progname_ = "calc_flux0";

    option long_options[] = {
        {"debug",       required_argument, NULL, 'd'},
        {"help",        required_argument, NULL, 'h'},
        {"verbose",     required_argument, NULL, 'v'},
        {0, 0, 0, 0}
    };

    // long option default

    SetOption(argc, argv, long_options);
    
    printf("ArgVal::Init: # of arg = %d\n", argc - optind);
    int narg = 6;
    if (argc - optind != narg){
        printf("# of arguments must be %d.\n", narg);
        Usage(stdout);
    }
    int iarg = optind;
    data_on_file_         = argv[iarg]; iarg++;
    phase_ratio_on_       = atof(argv[iarg]); iarg++;
    live_time_ratio_on_   = atof(argv[iarg]); iarg++;
    data_off_file_        = argv[iarg]; iarg++;
    phase_ratio_off_      = atof(argv[iarg]); iarg++;
    live_time_ratio_off_  = atof(argv[iarg]); iarg++;
}

void ArgValCalcFlux0::Print(FILE* fp) const
{
    fprintf(fp, "%s: g_flag_debug   : %d\n", __func__,
            g_flag_debug);
    fprintf(fp, "%s: g_flag_help    : %d\n", __func__,
            g_flag_help);
    fprintf(fp, "%s: g_flag_verbose : %d\n", __func__,
            g_flag_verbose);

    fprintf(fp, "%s: progname_       : %s\n", __func__,
            progname_.c_str());
    fprintf(fp, "%s: data_on_file_   : %s\n", __func__,
            data_on_file_.c_str());
    fprintf(fp, "%s: phase_ratio_on_ : %e\n", __func__,
            phase_ratio_on_);
    fprintf(fp, "%s: live_time_ratio_on_ : %e\n", __func__,
            live_time_ratio_on_);    
    fprintf(fp, "%s: data_off_file_  : %s\n", __func__,
            data_off_file_.c_str());
    fprintf(fp, "%s: phase_ratio_off_ : %e\n", __func__,
            phase_ratio_off_);
    fprintf(fp, "%s: live_time_ratio_off_ : %e\n", __func__,
            live_time_ratio_off_);
}


// private

void ArgValCalcFlux0::Null()
{
    progname_ = "";
    data_on_file_ = "";
    phase_ratio_on_ = 0.0;
    live_time_ratio_on_ = 0.0;
    data_off_file_ = "";
    phase_ratio_off_ = 0.0;
    live_time_ratio_off_ = 0.0;
}

void ArgValCalcFlux0::SetOption(int argc, char* argv[], option* long_options)
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


void ArgValCalcFlux0::Usage(FILE* fp) const
{
    fprintf(fp,
            "usage: %s [--help (0)] [--verbose (0)] [--debug (0)] "
            "data_on_file  phase_ratio_on  live_time_ratio_on  "
            "data_off_file  phase_ratio_off  live_time_ratio_off\n",
            progname_.c_str());
    fprintf(fp, "\n");
    fprintf(fp,
            "flux0 = flux_on - flux_off,\n");
    fprintf(fp,
            "flux_on = nevt_on / "
            "(phase_ratio_on * live_time_ratio_on)\n");
    fprintf(fp,
            "flux_off = nevt_off / "
            "(phase_ratio_off * live_time_ratio_off)\n");
    fprintf(fp, "if (data_off_file = none): "
            "flux0 = flux_on\n");
    abort();
}
