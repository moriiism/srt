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
data_list      = sys.argv[iarg]; iarg += 1
nfold          = int(sys.argv[iarg]); iarg += 1
outdir         = sys.argv[iarg]; iarg += 1
mu_list        = sys.argv[iarg]; iarg += 1
gamma_list     = sys.argv[iarg]; iarg += 1

print("data_list = ", data_list)
print("nfold = ", nfold)
print("outdir = ", outdir)
print("mu_list = ", mu_list)
print("gamma_list = ", gamma_list)

srt_dir = "/home/morii/work/github/moriiism/srt"

# data_list
print("data_list ...")
data_file_lst = []
phase_id_lst = []
phase_tag_lst = []
phase_ratio_lst = []
flux_0_lst = []
data_list_fptr = open(data_list, "r")
for line in data_list_fptr:
    line = line.rstrip()
    if(line[0] == "#"):
        continue
    print(line)
    (data_file, phase_id, phase_tag,
     phase_ratio, flux_0) = line.split()
    data_file_lst.append(data_file)
    phase_id_lst.append(phase_id)    
    phase_tag_lst.append(phase_tag)
    phase_ratio_lst.append(phase_ratio)
    flux_0_lst.append(flux_0)
    
data_list_fptr.close()
print(data_file_lst)
print(phase_id_lst)
print(phase_tag_lst)
print(phase_ratio_lst)
print(flux_0_lst)
print("data_list ... done.")

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

for gamma in gamma_par_lst:
    print("  gamma = ", gamma)
    gamma_dir = f"gamma{gamma:.1e}"
    gamma = float(gamma)

    print(" ", file=mu_rmse_file_fptr)
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

print("log x on", file=mu_rmse_file_fptr)
print("lw 7", file=mu_rmse_file_fptr)
print("line on", file=mu_rmse_file_fptr)
print("ma 6 on", file=mu_rmse_file_fptr)

mu_rmse_file_fptr.close()


#
#if __name__ == "__main__":
#    main()
#
