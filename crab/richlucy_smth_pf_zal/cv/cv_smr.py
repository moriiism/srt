#!/usr/bin/env python3
#
# richlucy_crab_cv_smr.py
#
'''
Overview:
  run after richlucy_crab_cv.py
Input:

Output:

Data:

Details:

'''
import os
import sys
import subprocess
import pandas as pd
import math

iarg = 1
rand_seed      = int(sys.argv[iarg]); iarg += 1
nfold          = int(sys.argv[iarg]); iarg += 1
outdir         = sys.argv[iarg]; iarg += 1
outfile_head   = sys.argv[iarg]; iarg += 1
mu_list_file   = sys.argv[iarg]; iarg += 1

print("rand_seed = ", rand_seed)
print("nfold = ", nfold)
print("outdir = ", outdir)
print("outfile_head = ", outfile_head)
print("mu_list_file = ", mu_list_file)


# mu_list
mu_lst = []
mu_list_file_fptr = open(mu_list_file, "r")
for line in mu_list_file_fptr:
    mu = line.rstrip('\n')
    mu_lst.append(mu)
mu_list_file_fptr.close()

print(mu_lst)

cmd = ["mkdir", outdir + "/" + "smr"]
print(cmd)
subprocess.call(cmd)
mu_rmse_file = f"{outdir}/smr/mu_rmse.qdp"
mu_rmse_file_fptr = open(mu_rmse_file, "w")
print("skip sing", file=mu_rmse_file_fptr)
print("read serr 2", file=mu_rmse_file_fptr)
print("! mu  rmse_ave  rmse_stddev", file=mu_rmse_file_fptr)

for mu in mu_lst:
    mu = float(mu)
    # eval
    for ifold in range(nfold):
        print("ifold = ", ifold)
        resp_file = f"{outdir}/resp/hxt_resp.fits"
        recfile = (f"{outdir}/rl_crab/mu{mu:.1e}/ifold{ifold:02}/{outfile_head}_pulsar+nebula.fits")
        valfile = (f"{outdir}/obs_cv/{outfile_head}_obs_{rand_seed:04}_" +
                   f"{nfold:02}fold{ifold:02}_vl.fits")
        outdir_eval = (f"{outdir}/rl_crab/mu{mu:.1e}/ifold{ifold:02}")
        outfile_head_eval = "eval"
        
        cmd = ["/home/morii/work/github/moriiism/srt/richlucy_crab/eval_val", 
               resp_file, recfile, valfile, str(nfold), outdir_eval, outfile_head_eval]
        print(cmd)
        subprocess.call(cmd)

        
    # average rmse
    rmse_sum = 0.0
    rmse_sum2 = 0.0
    nfold_exist = 0
    for ifold in range(nfold):
        rmse = 0.0
        rmse_file = f"{outdir}/rl_crab/mu{mu:.1e}/ifold{ifold:02}/eval_rmse.txt"
        if os.path.isfile(rmse_file) == False:
            print(f"warning: {rmse_file} does not exist.")
            continue
        else:
            nfold_exist += 1
        rmse_file_fptr = open(rmse_file, "r")
        for line in rmse_file_fptr:
            rmse = line.rstrip('\n')
        rmse_file_fptr.close()
        rmse_sum += float(rmse)
        rmse_sum2 += float(rmse) * float(rmse)

    if nfold_exist == 0:
        print(f"warning: nfold_exist == 0")
    else:
        if nfold_exist != nfold:
            print(f"warning: nfold_exist = {nfold_exist}, nfold = {nfold}")
        rmse_ave = rmse_sum / nfold_exist
        rmse_var = (rmse_sum2 - rmse_sum * rmse_sum / nfold_exist) / nfold_exist
        rmse_stddev = math.sqrt(rmse_var)
        print(f"{mu:.1e} {rmse_ave} {rmse_stddev}")
        print(f"{mu:.1e} {rmse_ave} {rmse_stddev}", file=mu_rmse_file_fptr)


print("log x on", file=mu_rmse_file_fptr)
print("lw 7", file=mu_rmse_file_fptr)
print("line on", file=mu_rmse_file_fptr)
print("ma 6 on", file=mu_rmse_file_fptr)

mu_rmse_file_fptr.close()


#
#if __name__ == "__main__":
#    main()

