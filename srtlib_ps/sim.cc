#include "sim.h"

// generate events from a probability distribution
void SrtlibSim::GenRandomEvtFromProbDist(const double* const prob_arr, int nbin,
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

void SrtlibSim::GenCVImageByPartition(int* const evt_arr, int nevt,
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
void SrtlibSim::GenCVImage(const double* const prob_arr, int nbin,
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
