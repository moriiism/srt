#include "arg_cvmap.h"

// public

void ArgValCvmap::Init(int argc, char* argv[])
{
    progname_ = "cvmap";

    option long_options[] = {
        {"debug",      required_argument, NULL, 'd'},
        {"help",       required_argument, NULL, 'h'},
        {"verbose",    required_argument, NULL, 'v'},
        {0, 0, 0, 0}
    };

    // long option default
    SetOption(argc, argv, long_options);

    if(0 < g_flag_verbose){
        printf("ArgVal::Init: # of arg = %d\n", argc - optind);
    }
    int narg = 6;
    if (argc - optind != narg){
        printf("# of arguments must be %d.\n", narg);
        Usage(stdout);
    }
    int iarg = optind;
    infile_       = argv[iarg]; iarg++;
    hist_info_file_  = argv[iarg]; iarg++;
    zrange_lo_    = atof(argv[iarg]); iarg++;
    zrange_up_    = atof(argv[iarg]); iarg++; 
    outdir_       = argv[iarg]; iarg++;
    outfile_head_ = argv[iarg]; iarg++;    
}

void ArgValCvmap::Print(FILE* fp) const
{
    fprintf(fp, "%s: g_flag_debug   : %d\n", __func__, g_flag_debug);
    fprintf(fp, "%s: g_flag_help    : %d\n", __func__, g_flag_help);
    fprintf(fp, "%s: g_flag_verbose : %d\n", __func__, g_flag_verbose);
    
    fprintf(fp, "%s: progname_      : %s\n", __func__, progname_.c_str());
    fprintf(fp, "%s: infile_        : %s\n", __func__, infile_.c_str());
    fprintf(fp, "%s: hist_info_file_ : %s\n", __func__, hist_info_file_.c_str());
    fprintf(fp, "%s: zrange_lo_     : %e\n", __func__, zrange_lo_);
    fprintf(fp, "%s: zrange_up_     : %e\n", __func__, zrange_up_);    
    fprintf(fp, "%s: outdir_        : %s\n", __func__, outdir_.c_str());
    fprintf(fp, "%s: outfile_head_  : %s\n", __func__, outfile_head_.c_str());    
}

void ArgValCvmap::Null()
{
    progname_     = "";
    infile_       = "";
    hist_info_file_ = "";
    zrange_lo_    = 0.0;
    zrange_up_    = 0.0;
    outdir_       = "";
    outfile_head_ = "";    
}

void ArgValCvmap::SetOption(int argc, char* argv[], option* long_options)
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
            printf("%s: error: getopt returned character code 0%o ??\n",
                   __func__, retopt);
            abort();
        }
    }
    if(0 < g_flag_verbose){
        MPrintInfo("done.");
    }
}

void ArgValCvmap::Usage(FILE* fp) const
{
    fprintf(fp,
            "usage: %s [--help (0)] [--verbose (0)] [--debug (0)] "
            "infile  hist_info_file  zrange_lo  zrange_up  outdir  outfile_head \n",
            progname_.c_str());
    exit(1);
}
