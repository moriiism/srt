#!/usr/bin/env python3
#
# cv_smr.py
#
'''
Overview:
  run after cv.py
Input:

Output:

Data:

Details:

'''
import os
import sys
import subprocess
import math

iarg = 1
nfold          = int(sys.argv[iarg]); iarg += 1
outdir         = sys.argv[iarg]; iarg += 1
mu_list        = sys.argv[iarg]; iarg += 1
gamma_list     = sys.argv[iarg]; iarg += 1

print("nfold = ", nfold)
print("outdir = ", outdir)
print("mu_list = ", mu_list)
print("gamma_list = ", gamma_list)

srt_dir = "/home/morii/work/github/moriiism/srt"

mu_par_lst = []
mu_list_fptr = open(mu_list, "r")
for line in mu_list_fptr:
    line = line.rstrip()
    if(line[0] == "#"):
        continue
    print(line)
    mu = line.rstrip('\n')
    mu = float(mu)
    mu_par_lst.append(mu)
mu_list_fptr.close()
print(mu_par_lst)

gamma_par_lst = []
gamma_list_fptr = open(gamma_list, "r")
for line in gamma_list_fptr:
    line = line.rstrip()
    if(line[0] == "#"):
        continue
    print(line)
    gamma = line.rstrip('\n')
    gamma = float(gamma)
    gamma_par_lst.append(gamma)
gamma_list_fptr.close()
print(gamma_par_lst)

cmd = ["mkdir", outdir + "/" + "smr"]
print(cmd)
subprocess.call(cmd)

mu_rmse_file = f"{outdir}/smr/mu_rmse.qdp"
mu_rmse_file_fptr = open(mu_rmse_file, "w")
print("skip sing", file=mu_rmse_file_fptr)
print("read serr 2", file=mu_rmse_file_fptr)
print("! mu  rmse_ave  rmse_stddev", file=mu_rmse_file_fptr)

rmse_ave_min = 1.0e+10 # very large
mu_min = -1.0
gamma_min = -1.0
rmse_stddev_min = -1.0
rmse_ave_plus_err_min = -1.0

flag_first = 1
for gamma in gamma_par_lst:
    print("  gamma = ", gamma)
    gamma_dir = f"gamma{gamma:.1e}"
    gamma = float(gamma)

    print(" ", file=mu_rmse_file_fptr)
    if flag_first == 1:
        flag_first = 0
    else:
        print("no", file=mu_rmse_file_fptr)        
    print(f"! gamma = {gamma:.1e}", gamma,
          file=mu_rmse_file_fptr)
    print(" ", file=mu_rmse_file_fptr)    

    for mu in mu_par_lst:
        print("mu = ", mu)
        mu_dir = f"mu{mu:.1e}"
        mu = float(mu)

        # average rmse
        rmse_sum = 0.0
        rmse_sum2 = 0.0
        for ifold in range(nfold):
            rmse_this = 0.0
            rmse_file = (
                f"{outdir}/rec/{mu_dir}/{gamma_dir}" + "/" +
                f"nfold{nfold:02}/ifold{ifold:02}" + "/" +
                "rl_rmse.txt")
            rmse_file_fptr = open(rmse_file, "r")
            for line in rmse_file_fptr:
                rmse_this = line.rstrip('\n')
            rmse_file_fptr.close()

            rmse_sum += float(rmse_this)
            rmse_sum2 += float(rmse_this) * float(rmse_this)
                
        rmse_ave = rmse_sum / nfold
        rmse_var = (rmse_sum2 - rmse_sum*rmse_sum/nfold)/nfold
        rmse_stddev = math.sqrt(rmse_var)
        print(f"{mu:.1e} {rmse_ave} {rmse_stddev} ! {gamma:.1e}")
        print(f"{mu:.1e} {rmse_ave} {rmse_stddev} ! {gamma:.1e}",
              file=mu_rmse_file_fptr)
        if rmse_ave < rmse_ave_min:
            rmse_ave_min = rmse_ave
            rmse_stddev_min = rmse_stddev
            mu_min = mu
            gamma_min = gamma
            rmse_ave_plus_err_min = rmse_ave + rmse_stddev 

print(f"rmse_ave_min = {rmse_ave_min} "
      f"at (mu, gamma) = ({mu_min}, {gamma_min})")
print(f"! rmse_ave_min = {rmse_ave_min} "
      f"at (mu, gamma) = ({mu_min}, {gamma_min})",
      file=mu_rmse_file_fptr)

print("la file", file=mu_rmse_file_fptr)
print("time off", file=mu_rmse_file_fptr)
print("la rot", file=mu_rmse_file_fptr)
print(f"la t \"{nfold}-fold Cross-Validation\"")
print("log x on", file=mu_rmse_file_fptr)
print("lw 7", file=mu_rmse_file_fptr)
print("line on", file=mu_rmse_file_fptr)
print("ma 6 on", file=mu_rmse_file_fptr)
print("csize 1.5", file=mu_rmse_file_fptr)
print("la x mu", file=mu_rmse_file_fptr)
print("la y RMSE", file=mu_rmse_file_fptr)

print(f"! rmse_ave_min({rmse_ave_min}) " +
      f" rmse_stddev_min({rmse_stddev_min}) = " +
      f"{rmse_ave_plus_err_min}",
      file=mu_rmse_file_fptr)
print(f"label 1 pos {mu_min} {rmse_ave_plus_err_min} " +
      "line 0 1.0 \" \" lst 4",
      file=mu_rmse_file_fptr)
print("hard mu_rmse.ps/cps", file=mu_rmse_file_fptr)

mu_rmse_file_fptr.close()


#
#if __name__ == "__main__":
#    main()
#
