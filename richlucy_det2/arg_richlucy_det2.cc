#include "arg_richlucy_det2.h"

// public

void ArgValRichlucyDet2::Init(int argc, char* argv[])
{
    progname_ = "richlucy_det2";

    option long_options[] = {
        {"debug",       required_argument, NULL, 'd'},
        {"help",        required_argument, NULL, 'h'},
        {"verbose",     required_argument, NULL, 'v'},
        {0, 0, 0, 0}
    };

    // long option default

    SetOption(argc, argv, long_options);
    
    printf("ArgVal::Init: # of arg = %d\n", argc - optind);
    int narg = 16;
    if (argc - optind != narg){
        printf("# of arguments must be %d.\n", narg);
        Usage(stdout);
    }
    int iarg = optind;
    datafile1_       = argv[iarg]; iarg++;
    datafile2_       = argv[iarg]; iarg++;    
    bg_file1_        = argv[iarg]; iarg++;
    bg_file2_        = argv[iarg]; iarg++;    
    resp_norm_file1_ = argv[iarg]; iarg++;
    resp_norm_file2_ = argv[iarg]; iarg++;    
    eff_file_        = argv[iarg]; iarg++;
    nskyx_           = atoi(argv[iarg]); iarg++;
    nskyy_           = atoi(argv[iarg]); iarg++;
    ndetx_           = atoi(argv[iarg]); iarg++;
    ndety_           = atoi(argv[iarg]); iarg++;
    outdir_          = argv[iarg]; iarg++;
    outfile_head_    = argv[iarg]; iarg++;
    nem_             = atoi(argv[iarg]); iarg++;
    tol_em_          = atof(argv[iarg]); iarg++;
    acc_method_      = argv[iarg]; iarg++;
}

void ArgValRichlucyDet2::Print(FILE* fp) const
{
    fprintf(fp, "%s: g_flag_debug   : %d\n",
            __func__, g_flag_debug);
    fprintf(fp, "%s: g_flag_help    : %d\n",
            __func__, g_flag_help);
    fprintf(fp, "%s: g_flag_verbose : %d\n",
            __func__, g_flag_verbose);

    fprintf(fp, "%s: progname_       : %s\n",
            __func__, progname_.c_str());
    fprintf(fp, "%s: datafile1_       : %s\n",
            __func__, datafile1_.c_str());
    fprintf(fp, "%s: datafile2_       : %s\n",
            __func__, datafile2_.c_str());    
    fprintf(fp, "%s: bg_file1_       : %s\n",
            __func__, bg_file1_.c_str());
    fprintf(fp, "%s: bg_file2_       : %s\n",
            __func__, bg_file2_.c_str());        
    fprintf(fp, "%s: resp_norm_file1_      : %s\n",
            __func__, resp_norm_file1_.c_str());
    fprintf(fp, "%s: resp_norm_file2_      : %s\n",
            __func__, resp_norm_file2_.c_str());
    fprintf(fp, "%s: eff_file_       : %s\n",
            __func__, eff_file_.c_str());
    fprintf(fp, "%s: nskyx_          : %d\n",
            __func__, nskyx_);
    fprintf(fp, "%s: nskyy_          : %d\n",
            __func__, nskyy_);
    fprintf(fp, "%s: ndetx_          : %d\n",
            __func__, ndetx_);
    fprintf(fp, "%s: ndety_          : %d\n",
            __func__, ndety_);
    fprintf(fp, "%s: outdir_         : %s\n",
            __func__, outdir_.c_str());
    fprintf(fp, "%s: outfile_head_   : %s\n",
            __func__, outfile_head_.c_str());
    fprintf(fp, "%s: nem_            : %d\n",
            __func__, nem_);
    fprintf(fp, "%s: tol_em_         : %f\n",
            __func__, tol_em_);
    fprintf(fp, "%s: acc_method_     : %s\n",
            __func__, acc_method_.c_str());
}

// private

void ArgValRichlucyDet2::Null()
{
    progname_ = "";
    datafile1_ = "";
    datafile2_ = "";    
    bg_file1_ = "";
    bg_file2_ = "";    
    resp_norm_file1_  = "";
    resp_norm_file2_  = "";    
    eff_file_  = "";
    nskyx_     = 0;
    nskyy_     = 0;
    ndetx_     = 0;
    ndety_     = 0;
    outdir_   = "";
    outfile_head_ = "";
    nem_      = 0;    
    tol_em_   = 0.0;
    acc_method_ = "";
}

void ArgValRichlucyDet2::SetOption(int argc, char* argv[], option* long_options)
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


void ArgValRichlucyDet2::Usage(FILE* fp) const
{
    fprintf(fp,
            "usage: %s [--help (0)] [--verbose (0)] [--debug (0)] "
            "datafile1  datafile2  bg_file1  bg_file2  "
            "resp_norm_file1  resp_norm_file2  eff_file  "
            "nskyx  nskyy  ndetx  ndety  "
            "outdir  outfile_head  "
            "nem  tol_em  acc_method\n",
            progname_.c_str());
    abort();
}
