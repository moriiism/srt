#include "arg_imgadd.h"

// public

void ArgValImgadd::Init(int argc, char* argv[])
{
    progname_ = "imgadd";

    option long_options[] = {
        {"debug",       required_argument, NULL, 'd'},
        {"help",        required_argument, NULL, 'h'},
        {"verbose",     required_argument, NULL, 'v'},
        {0, 0, 0, 0}
    };

    // long option default

    SetOption(argc, argv, long_options);
    
    printf("ArgVal::Init: # of arg = %d\n", argc - optind);
    int narg = 5;
    if (argc - optind != narg){
        printf("# of arguments must be %d.\n", narg);
        Usage(stdout);
    }
    int iarg = optind;
    imgfile1_       = argv[iarg]; iarg++;
    imgfile2_       = argv[iarg]; iarg++;
    coeff1_         = atof(argv[iarg]); iarg++;
    coeff2_         = atof(argv[iarg]); iarg++;
    outfile_        = argv[iarg]; iarg++;
}

void ArgValImgadd::Print(FILE* fp) const
{
    fprintf(fp, "%s: g_flag_debug   : %d\n", __func__, g_flag_debug);
    fprintf(fp, "%s: g_flag_help    : %d\n", __func__, g_flag_help);
    fprintf(fp, "%s: g_flag_verbose : %d\n", __func__, g_flag_verbose);

    fprintf(fp, "%s: progname_       : %s\n", __func__, progname_.c_str());
    fprintf(fp, "%s: imgfile1_       : %s\n", __func__, imgfile1_.c_str());
    fprintf(fp, "%s: imgfile2_       : %s\n", __func__, imgfile2_.c_str());
    fprintf(fp, "%s: coeff1_         : %e\n", __func__, coeff1_);
    fprintf(fp, "%s: coeff2_         : %e\n", __func__, coeff2_);
    fprintf(fp, "%s: outfile_        : %s\n", __func__, outfile_.c_str());
}

// private

void ArgValImgadd::Null()
{
    progname_  = "";
    imgfile1_  = "";
    imgfile2_  = "";
    coeff1_    = 0.0;
    coeff2_    = 0.0;    
    outfile_   = "";
}

void ArgValImgadd::SetOption(int argc, char* argv[], option* long_options)
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


void ArgValImgadd::Usage(FILE* fp) const
{
    fprintf(fp,
            "usage: %s [--help (0)] [--verbose (0)] [--debug (0)] "
            "imgfile1  imgfile2  coeff1  coeff2  outfile\n",
            progname_.c_str());
    abort();
}
