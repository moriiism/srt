#include "mir_math.h"
#include "mif_fits.h"
#include "mif_img_info.h"
#include "mi_time.h"
#include "arg_imgadd.h"

// global variable 
int g_flag_debug = 0;
int g_flag_help = 0;
int g_flag_verbose = 0;

int main(int argc, char* argv[])
{
    int status_prog = kRetNormal;
    
    ArgValImgadd* argval = new ArgValImgadd;
    argval->Init(argc, argv);
    argval->Print(stdout);

    int nskyx1 = MifFits::GetAxisSize(argval->GetImgfile1(), 0);
    int nskyy1 = MifFits::GetAxisSize(argval->GetImgfile1(), 1);    
    int nskyx2 = MifFits::GetAxisSize(argval->GetImgfile2(), 0);
    int nskyy2 = MifFits::GetAxisSize(argval->GetImgfile2(), 1);
    if( (nskyx1 != nskyx2) ||
        (nskyy1 != nskyy2) ){
        printf("nskyx1 != nskyx2 or nskyy1 != nskyy2\n");
        abort();
    }

    // load img1
    double* img1_arr = NULL;
    MifImgInfo* img_info_1 = new MifImgInfo;
    img_info_1->InitSetImg(1, 1, nskyx1, nskyy1);
    int bitpix_img1 = 0;
    MifFits::InFitsImageD(argval->GetImgfile1(), img_info_1,
                          &bitpix_img1, &img1_arr);
    delete img_info_1;

    // load img2
    double* img2_arr = NULL;
    MifImgInfo* img_info_2 = new MifImgInfo;
    img_info_2->InitSetImg(1, 1, nskyx2, nskyy2);
    int bitpix_img2 = 0;
    MifFits::InFitsImageD(argval->GetImgfile2(), img_info_2,
                          &bitpix_img2, &img2_arr);
    delete img_info_2;

    int nsky = nskyx1 * nskyy1;
    double* add_arr = new double[nsky];
    for(int isky = 0; isky < nsky; isky++){
        add_arr[isky] = img1_arr[isky] * argval->GetCoeff1()
            + img2_arr[isky] * argval->GetCoeff2();
    }


    long naxes[2];
    naxes[0] = nskyx1;
    naxes[1] = nskyy1;
    int bitpix = -64;
    MifFits::OutFitsImageD(argval->GetOutfile(),
                           2, bitpix,
                           naxes, add_arr);
    
    return status_prog;
}
