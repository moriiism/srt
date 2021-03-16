#include "arg_richlucy_bg.h"

// public

void ArgValRichlucyBg::Init(int argc, char* argv[])
{
    progname_ = "richlucy";

    option long_options[] = {
        {"debug",       required_argument, NULL, 'd'},
        {"help",        required_argument, NULL, 'h'},
        {"verbose",     required_argument, NULL, 'v'},
        {0, 0, 0, 0}
    };

    // long option default

    SetOption(argc, argv, long_options);
    
    printf("ArgVal::Init: # of arg = %d\n", argc - optind);
    int narg = 13;
    if (argc - optind != narg){
        printf("# of arguments must be %d.\n", narg);
        Usage(stdout);
    }
    int iarg = optind;
    respdir_        = argv[iarg]; iarg++;
    datafile_       = argv[iarg]; iarg++;
    skyfile_        = argv[iarg]; iarg++;
    bgfile_         = argv[iarg]; iarg++;
    outdir_         = argv[iarg]; iarg++;
    outfile_head_   = argv[iarg]; iarg++;
    nloop_main_     = atoi(argv[iarg]); iarg++;
    nloop_em_       = atoi(argv[iarg]); iarg++;
    nloop_newton_   = atoi(argv[iarg]); iarg++;
    tol_main_       = atof(argv[iarg]); iarg++;
    tol_em_         = atof(argv[iarg]); iarg++;
    tol_newton_     = atof(argv[iarg]); iarg++;
    epsilon_        = atof(argv[iarg]); iarg++;
}

void ArgValRichlucyBg::Print(FILE* fp) const
{
    fprintf(fp, "%s: g_flag_debug   : %d\n", __func__, g_flag_debug);
    fprintf(fp, "%s: g_flag_help    : %d\n", __func__, g_flag_help);
    fprintf(fp, "%s: g_flag_verbose : %d\n", __func__, g_flag_verbose);

    fprintf(fp, "%s: progname_       : %s\n", __func__, progname_.c_str());
    fprintf(fp, "%s: respdir_        : %s\n", __func__, respdir_.c_str());
    fprintf(fp, "%s: datafile_       : %s\n", __func__, datafile_.c_str());
    fprintf(fp, "%s: skyfile_        : %s\n", __func__, skyfile_.c_str());
    fprintf(fp, "%s: bgfile_         : %s\n", __func__, bgfile_.c_str());
    fprintf(fp, "%s: outdir_         : %s\n", __func__, outdir_.c_str());
    fprintf(fp, "%s: outfile_head_   : %s\n", __func__, outfile_head_.c_str());
    fprintf(fp, "%s: nloop_main_     : %d\n", __func__, nloop_main_);
    fprintf(fp, "%s: nloop_em_       : %d\n", __func__, nloop_em_);
    fprintf(fp, "%s: nloop_newton_   : %d\n", __func__, nloop_newton_);
    fprintf(fp, "%s: tol_main_       : %f\n", __func__, tol_main_);
    fprintf(fp, "%s: tol_em_         : %f\n", __func__, tol_em_);
    fprintf(fp, "%s: tol_newton_     : %f\n", __func__, tol_newton_);
    fprintf(fp, "%s: epsilon_        : %f\n", __func__, epsilon_);
}

// private

void ArgValRichlucyBg::Null()
{
    progname_ = "";
    respdir_  = "";
    datafile_ = "";
    skyfile_  = "";
    bgfile_   = "";
    outdir_   = "";
    outfile_head_ = "";
    nloop_main_   = 0;    
    nloop_em_     = 0;
    nloop_newton_ = 0;
    tol_main_     = 0.0;
    tol_em_       = 0.0;
    tol_newton_   = 0.0;
    epsilon_      = 0.0;
}

void ArgValRichlucyBg::SetOption(int argc, char* argv[], option* long_options)
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


void ArgValRichlucyBg::Usage(FILE* fp) const
{
    fprintf(fp,
            "usage: %s [--help (0)] [--verbose (0)] [--debug (0)] "
            "respdir  datafile  skyfile  bgfile  "
            "outdir  outfile_head  nloop_main  nloop_em  nloop_newton "
            "tol_main  tol_em  tol_newton  epsilon\n",
            progname_.c_str());
    abort();
}
