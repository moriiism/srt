#include "arg_npeak.h"

// public

void ArgValNpeak::Init(int argc, char* argv[])
{
    progname_ = "npeak";

    option long_options[] = {
        {"debug",       required_argument, NULL, 'd'},
        {"help",        required_argument, NULL, 'h'},
        {"verbose",     required_argument, NULL, 'v'},
        {0, 0, 0, 0}
    };

    // long option default

    SetOption(argc, argv, long_options);
    
    printf("ArgVal::Init: # of arg = %d\n", argc - optind);
    int narg = 12;
    if (argc - optind != narg){
        printf("# of arguments must be %d.\n", narg);
        Usage(stdout);
    }
    int iarg = optind;
    infile_         = argv[iarg]; iarg++;
    xpos_           = atof(argv[iarg]); iarg++;
    ypos_           = atof(argv[iarg]); iarg++;
    theta_          = atof(argv[iarg]); iarg++;
    y_lo_           = atof(argv[iarg]); iarg++;
    y_up_           = atof(argv[iarg]); iarg++;
    nbinx_new_      = atoi(argv[iarg]); iarg++;
    xlo_new_        = atof(argv[iarg]); iarg++;
    xup_new_        = atof(argv[iarg]); iarg++;
    significance_   = atof(argv[iarg]); iarg++;
    outdir_         = argv[iarg]; iarg++;
    outfile_head_   = argv[iarg]; iarg++;
}

void ArgValNpeak::Print(FILE* fp) const
{
    fprintf(fp, "%s: g_flag_debug   : %d\n", __func__, g_flag_debug);
    fprintf(fp, "%s: g_flag_help    : %d\n", __func__, g_flag_help);
    fprintf(fp, "%s: g_flag_verbose : %d\n", __func__, g_flag_verbose);

    fprintf(fp, "%s: progname_       : %s\n", __func__, progname_.c_str());
    fprintf(fp, "%s: infile_         : %s\n", __func__, infile_.c_str());
    fprintf(fp, "%s: xpos_           : %e\n", __func__, xpos_);
    fprintf(fp, "%s: ypos_           : %e\n", __func__, ypos_);
    fprintf(fp, "%s: theta_          : %e\n", __func__, theta_);
    fprintf(fp, "%s: y_lo_           : %e\n", __func__, y_lo_);
    fprintf(fp, "%s: y_up_           : %e\n", __func__, y_up_);
    fprintf(fp, "%s: nbinx_new_      : %d\n", __func__, nbinx_new_);
    fprintf(fp, "%s: xlo_new_        : %e\n", __func__, xlo_new_);
    fprintf(fp, "%s: xup_new_        : %e\n", __func__, xup_new_);
    fprintf(fp, "%s: significance_   : %e\n", __func__, significance_);
    fprintf(fp, "%s: outdir_         : %s\n", __func__, outdir_.c_str());
    fprintf(fp, "%s: outfile_head_   : %s\n", __func__, outfile_head_.c_str());
}

// private

void ArgValNpeak::Null()
{
    progname_  = "";
    infile_    = "";
    xpos_      = 0.0;
    ypos_      = 0.0;
    theta_     = 0.0;
    y_lo_      = 0.0;
    y_up_      = 0.0;
    nbinx_new_ = 0;
    xlo_new_   = 0.0;
    xup_new_   = 0.0;
    significance_ = 0.0;
    outdir_    = "";
    outfile_head_   = "";
}

void ArgValNpeak::SetOption(int argc, char* argv[], option* long_options)
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


void ArgValNpeak::Usage(FILE* fp) const
{
    fprintf(fp,
            "usage: %s [--help (0)] [--verbose (0)] [--debug (0)] "
            "infile  xpos  ypos  theta  y_lo  y_up  nbinx_new  xlo_new  xup_new  significance  outdir  outfile_head\n",
            progname_.c_str());
    abort();
}
