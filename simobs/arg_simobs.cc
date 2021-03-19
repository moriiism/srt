#include "arg_simobs.h"

// public

void ArgValSimobs::Init(int argc, char* argv[])
{
    progname_ = "simobs";

    option long_options[] = {
        {"debug",       required_argument, NULL, 'd'},
        {"help",        required_argument, NULL, 'h'},
        {"verbose",     required_argument, NULL, 'v'},
        {0, 0, 0, 0}
    };

    // long option default

    SetOption(argc, argv, long_options);
    
    printf("ArgVal::Init: # of arg = %d\n", argc - optind);
    int narg = 8;
    if (argc - optind != narg){
        printf("# of arguments must be %d.\n", narg);
        Usage(stdout);
    }
    int iarg = optind;
    respfile_       = argv[iarg]; iarg++;
    srcfile_        = argv[iarg]; iarg++;
    nevt_src_       = atoi(argv[iarg]); iarg++;
    bgfile_         = argv[iarg]; iarg++;
    nevt_bg_        = atoi(argv[iarg]); iarg++;
    rand_seed_det_  = atoi(argv[iarg]); iarg++;
    outdir_         = argv[iarg]; iarg++;
    outfile_head_   = argv[iarg]; iarg++;
}

void ArgValSimobs::Print(FILE* fp) const
{
    fprintf(fp, "%s: g_flag_debug   : %d\n", __func__, g_flag_debug);
    fprintf(fp, "%s: g_flag_help    : %d\n", __func__, g_flag_help);
    fprintf(fp, "%s: g_flag_verbose : %d\n", __func__, g_flag_verbose);

    fprintf(fp, "%s: progname_       : %s\n", __func__, progname_.c_str());
    fprintf(fp, "%s: respfile_       : %s\n", __func__, respfile_.c_str());
    fprintf(fp, "%s: srcfile_        : %s\n", __func__, srcfile_.c_str());
    fprintf(fp, "%s: nevt_src_       : %d\n", __func__, nevt_src_);
    fprintf(fp, "%s: bgfile_         : %s\n", __func__, bgfile_.c_str());
    fprintf(fp, "%s: nevt_bg_        : %d\n", __func__, nevt_bg_);
    fprintf(fp, "%s: rand_seed_det_  : %d\n", __func__, rand_seed_det_);
    fprintf(fp, "%s: outdir_         : %s\n", __func__, outdir_.c_str());
    fprintf(fp, "%s: outfile_head_   : %s\n", __func__, outfile_head_.c_str());
}

// private

void ArgValSimobs::Null()
{
    progname_  = "";
    respfile_  = "";
    srcfile_   = "";
    nevt_src_  = 0;
    bgfile_    = "";
    nevt_bg_   = 0;
    rand_seed_det_ = 0;
    outdir_    = "";
    outfile_head_ = "";
}

void ArgValSimobs::SetOption(int argc, char* argv[], option* long_options)
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


void ArgValSimobs::Usage(FILE* fp) const
{
    fprintf(fp,
            "usage: %s [--help (0)] [--verbose (0)] [--debug (0)] "
            "respfile  srcfile  nevt_src  bgfile  nevt_bg  "
            "rand_seed_det  outdir  outfile_head\n",
            progname_.c_str());
    abort();
}
