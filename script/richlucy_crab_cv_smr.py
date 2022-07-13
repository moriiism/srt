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
nfold          = int(sys.argv[iarg]); iarg += 1
outdir         = sys.argv[iarg]; iarg += 1
mu_list_file   = sys.argv[iarg]; iarg += 1

print("nfold = ", nfold)
print("outdir = ", outdir)
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
mu_helldist_file = f"{outdir}/smr/mu_helldist.qdp"
mu_helldist_file_fptr = open(mu_helldist_file, "w")
print("skip sing")
print("read serr 2")
print("! mu  helldist_ave  helldist_stddev",
      file=mu_helldist_file_fptr)

for mu in mu_lst:
    mu = float(mu)
    # average helldist
    helldist_sum = 0.0
    helldist_sum2 = 0.0
    nfold_exist = 0
    for ifold in range(nfold):
        helldist = 0.0
        helldist_file = (
            f"{outdir}/rl_crab/mu{mu:.1e}/ifold{ifold:02}/eval_helldist.txt")
        if os.path.isfile(helldist_file) == False:
            print(f"warning: {helldist_file} does not exist.")
            continue
        else:
            nfold_exist += 1
        helldist_file_fptr = open(helldist_file, "r")
        for line in helldist_file_fptr:
            helldist = line.rstrip('\n')
        helldist_file_fptr.close()
        helldist_sum += float(helldist)
        helldist_sum2 += float(helldist) * float(helldist)

    if nfold_exist == 0:
        print(f"warning: nfold_exist == 0")
    else:
        if nfold_exist != nfold:
            print(f"warning: nfold_exist = {nfold_exist}, nfold = {nfold}")
        helldist_ave = helldist_sum / nfold_exist
        helldist_var = (helldist_sum2 - helldist_sum * helldist_sum / nfold_exist) / nfold_exist
        helldist_stddev = math.sqrt(helldist_var)
        print(helldist_ave)
        print(helldist_stddev)
        print(f"{mu:.1e} {helldist_ave} {helldist_stddev}",
              file=mu_helldist_file_fptr)

mu_helldist_file_fptr.close()


#
#if __name__ == "__main__":
#    main()

