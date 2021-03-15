#include "sub.h"

void GetNdet(string respdir, int* const ndetx_ptr, int* const ndety_ptr)
{
    char infile[kLineSize];
    sprintf(infile, "%s/gimage_000_000.img", respdir.c_str());
    int naxis = MifFits::GetNaxis(infile);
    printf("GetNdet: naxis = %d\n", naxis);
    if(2 != naxis){
        printf("GetNdet: bad naxis = %d\n", naxis);
        abort();
    }
    int* ndet_arr = new int[naxis];
    for(int iaxis = 0; iaxis < naxis; iaxis ++){
        ndet_arr[iaxis] = MifFits::GetAxisSize(infile, iaxis);
    }
    *ndetx_ptr = ndet_arr[0];
    *ndety_ptr = ndet_arr[1];
    delete [] ndet_arr;
}

void LoadResp(string respdir, int nskyx, int nskyy,
              double** const mat_arr_ptr,
              int* const ndetx_ptr,
              int* const ndety_ptr)
{
    printf("--- LoadResp --- \n");

//    row: detector
//    col: sky
    
    int ndetx = 0;
    int ndety = 0;
    GetNdet(respdir, &ndetx, &ndety);
    printf("LoadResp: ndetx = %d\n", ndetx);
    printf("LoadResp: ndety = %d\n", ndety);
    MifImgInfo* img_info = new MifImgInfo;
    img_info->InitSetImg(1, 1, ndetx, ndety);

    int ndet = ndetx * ndety;
    int nsky = nskyx * nskyy;
    printf("LoadResp: ndet = %d, nsky = %d\n", ndet, nsky);
    double* mat_arr = new double [ndet * nsky];
    for(int iskyy = 0; iskyy < nskyy; iskyy ++){
        for(int iskyx = 0; iskyx < nskyx; iskyx ++){
            char infile[kLineSize];
            sprintf(infile, "%s/gimage_%3.3d_%3.3d.img", respdir.c_str(), iskyx, iskyy);

            int bitpix = 0;
            double* data_arr = NULL;
            MifFits::InFitsImageD(infile, img_info,
                                  &bitpix, &data_arr);
            double sum = 0.0;
            for(int idet = 0; idet < ndet; idet ++){
                sum += data_arr[idet];
            }
            
            int isky = nskyx * iskyy + iskyx;
            int imat = isky * ndet;
            for(int idet = 0; idet < ndet; idet ++){
                mat_arr[imat + idet] = data_arr[idet] / sum;
            }
            delete [] data_arr;
        }
    }
    delete img_info;
    
    *mat_arr_ptr = mat_arr;
    *ndetx_ptr = ndetx;
    *ndety_ptr = ndety;
}

// generate events from a probability distribution
void GenRandomEvtFromProbDist(const double* const prob_arr, int nbin,
                              int nevt, int rand_seed,
                              double* const out_bin_arr,
                              int*   const out_evt_arr)
{
    for(int ibin = 0; ibin < nbin; ibin ++){
        out_bin_arr[ibin] = 0.0;
    }
    for(int ievt = 0; ievt < nevt; ievt ++){
        out_evt_arr[ievt] = 0;
    }
    
    // cumulative dist
    double* cum_arr = new double [nbin];
    double cum = 0.0;
    for(int ibin = 0; ibin < nbin; ibin ++){
        cum += prob_arr[ibin];
        cum_arr[ibin] = cum;
    }
    TRandom3* trand = new TRandom3(rand_seed);
    for(int ievt = 0; ievt < nevt; ievt++){
        double rand = trand->Rndm();
        int ibin_find = MiSort::BinarySearch(nbin, cum_arr, rand);
        out_bin_arr[ibin_find + 1] ++;
        out_evt_arr[ievt] = ibin_find + 1;
    }
    delete trand;
    delete [] cum_arr;
}

void GenCVImageByPartition(int* const evt_arr, int nevt,
                           int nfold, int rand_seed_partition,
                           int nbin,
                           double** const out_tr_arr,
                           double** const out_vl_arr)
{
    for(int ifold = 0; ifold < nfold; ifold ++){
        for(int ibin = 0; ibin < nbin; ibin ++){
            out_tr_arr[ifold][ibin] = 0.0;
            out_vl_arr[ifold][ibin] = 0.0;
        }
    }

    TRandom3* trand = new TRandom3(rand_seed_partition);
    double* rand_arr = new double[nevt];
    for(int ievt = 0; ievt < nevt; ievt ++){
        double rand = trand->Rndm();
        rand_arr[ievt] = rand;
    }
    int* index_arr = new int[nevt];
    MiSort::Sort(nevt, rand_arr, index_arr, 1);

    double nevt_fold = double(nevt) / double(nfold);
    for(int ievt = 0; ievt < nevt; ievt ++){
        int ifold = (int) floor(index_arr[ievt] / nevt_fold);
        out_vl_arr[ifold][ evt_arr[ievt] ] ++;
        for(int ifold_tr = 0; ifold_tr < nfold; ifold_tr ++){
            if(ifold_tr != ifold){
                out_tr_arr[ifold_tr][ evt_arr[ievt] ] ++;
            }
        }
    }
    delete trand;
    delete [] rand_arr;
    delete [] index_arr;
}


// Make images for N-fold cross-validation
void GenCVImage(const double* const prob_arr, int nbin,
                int nevt, int rand_seed, int nfold,
                double** const out_tr_arr,
                double** const out_vl_arr,
                double*  const out_arr)
{
    for(int ibin = 0; ibin < nbin; ibin ++){
        out_arr[ibin] = 0.0;
    }
    for(int ifold = 0; ifold < nfold; ifold ++){
        for(int ibin = 0; ibin < nbin; ibin ++){
            out_tr_arr[ifold][ibin] = 0.0;
            out_vl_arr[ifold][ibin] = 0.0;
        }
    }
    // cumulative dist
    double* cum_arr = new double [nbin];
    double cum = 0.0;
    for(int ibin = 0; ibin < nbin; ibin ++){
        cum += prob_arr[ibin];
        cum_arr[ibin] = cum;
    }

    int nevt_fold = nevt / nfold;
    TRandom3* trand = new TRandom3(rand_seed);
    for(int ievt = 0; ievt < nevt; ievt++){
        double rand = trand->Rndm();
        int ibin_find = MiSort::BinarySearch(nbin, cum_arr, rand);
        out_arr[ibin_find + 1] ++;
        int ifold = ievt / nevt_fold;
        out_vl_arr[ifold][ibin_find + 1] ++;
        for(int ifold_tr = 0; ifold_tr < nfold; ifold_tr ++){
            if(ifold_tr != ifold){
                out_tr_arr[ifold_tr][ibin_find + 1] ++;
            }
        }
    }
    delete trand;
    delete [] cum_arr;
}
