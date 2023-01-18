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
#import math

iarg = 1
data_list      = sys.argv[iarg]; iarg += 1
rand_seed      = int(sys.argv[iarg]); iarg += 1
nfold          = int(sys.argv[iarg]); iarg += 1
outdir         = sys.argv[iarg]; iarg += 1
outfile_head   = sys.argv[iarg]; iarg += 1
mu_list        = sys.argv[iarg]; iarg += 1
gamma_list     = sys.argv[iarg]; iarg += 1

print("data_list = ", data_list)
print("rand_seed = ", rand_seed)
print("nfold = ", nfold)
print("outdir = ", outdir)
print("outfile_head = ", outfile_head)
print("mu_list = ", mu_list)
print("gamma_list = ", gamma_list)

srt_dir = "/home/morii/work/github/moriiism/srt"

# data_list
print("data_list ...")
data_file_lst = []
phase_tag_lst = []
phase_ratio_lst = []
flux_0_lst = []
data_list_fptr = open(data_list, "r")
for line in data_list_fptr:
    line = line.rstrip()
    if(line[0] == "#"):
        continue
    print(line)
    (data_file, phase_tag, phase_ratio, flux_0) = line.split()
    data_file_lst.append(data_file)
    phase_tag_lst.append(phase_tag)
    phase_ratio_lst.append(phase_ratio)
    flux_0_lst.append(flux_0)
    
data_list_fptr.close()
print(data_file_lst)
print(phase_tag_lst)
print(phase_ratio_lst)
print(flux_0_lst)
print("data_list ... done.")

# mu_list
mu_lst = []
mu_list_fptr = open(mu_list, "r")
for line in mu_list_fptr:
    line = line.rstrip()
    if(line[0] == "#"):
        continue
    print(line)
    mu = float(line)
    mu_lst.append(mu)
mu_list_fptr.close()
print(mu_lst)

# gamma_list
gamma_lst = []
gamma_list_fptr = open(gamma_list, "r")
for line in gamma_list_fptr:
    line = line.rstrip()
    if(line[0] == "#"):
        continue
    print(line)
    gamma = float(line)
    gamma_lst.append(gamma)
gamma_list_fptr.close()
print(gamma_lst)


cmd = ["mkdir", outdir + "/" + "smr"]
print(cmd)
subprocess.call(cmd)
mu_rmse_file = f"{outdir}/smr/mu_rmse.qdp"
mu_rmse_file_fptr = open(mu_rmse_file, "w")
print("skip sing", file=mu_rmse_file_fptr)
print("read serr 2", file=mu_rmse_file_fptr)
print("! mu  rmse_ave  rmse_stddev", file=mu_rmse_file_fptr)

for mu in mu_lst:
    print("mu = %e", mu)
    mu_dir = f"mu{mu:.1e}"
    mu = float(mu)
    for gamma in gamma_lst:
        print("  gamma = %s", gamma)
        gamma_dir = f"gamma{gamma:.1e}"
        gamma = float(gamma)
        for ifold in range(nfold):
            print("ifold = ", ifold)

            rmse_sum2 = 0.0
            for index in range(len(phase_tag_lst)):
                # eval
                resp_file = f"{outdir}/resp/hxt1_resp_norm.fits"
                recfile = (
                    f"{outdir}/rec/{mu_dir}/{gamma_dir}" + "/" +
                    f"nfold{nfold:02}/ifold{ifold:02}" + "/" +
                    f"rec_rec_{index:02}.fits")
                valfile = (
                    f"{outdir}/obs_cv/{phase_tag_lst[index]}" +
                    "/" + f"{outfile_head}_obs_{rand_seed:04}_" +
                    f"{nfold:02}fold{ifold:02}_vl.fits")
                outdir_eval = (
                    f"{outdir}/rec/{mu_dir}/{gamma_dir}" + "/" +
                    f"nfold{nfold:02}/ifold{ifold:02}")
                outfile_head_eval = f"eval_{index:02}"

                cmd = [srt_dir + "/" + "crab/richlucy_smth_pf_zal"
                       + "/" + "eval_val",
                       resp_file, recfile, valfile, str(nfold),
                       outdir_eval, outfile_head_eval]
                print(cmd)
                subprocess.call(cmd)

                # read eval_rmse.txt
                rmse_this = 0.0
                rmse_file = (
                    f"{outdir}/rec/{mu_dir}/{gamma_dir}" + "/" +
                    f"nfold{nfold:02}/ifold{ifold:02}" + "/" +
                    f"eval_{index:02}_rmse.txt")
                rmse_file_fptr = open(rmse_file, "r")
                for line in rmse_file_fptr:
                    rmse_this = line.rstrip('\n')
                rmse_file_fptr.close()
                rmse_sum2 += float(rmse_this) * float(rmse_this)
                
            # add eval.txt
            rmse = math.sqrt(rmse_sum2)
            
##################################
exit(1)

#
#        # average rmse
#        rmse_sum = 0.0
#        rmse_sum2 = 0.0
#        nfold_exist = 0
#        for ifold in range(nfold):
#            rmse = 0.0
#        rmse_file = f"{outdir}/rl_crab/mu{mu:.1e}/ifold{ifold:02}/eval_rmse.txt"
#        if os.path.isfile(rmse_file) == False:
#            print(f"warning: {rmse_file} does not exist.")
#            continue
#        else:
#            nfold_exist += 1
#        rmse_file_fptr = open(rmse_file, "r")
#        for line in rmse_file_fptr:
#            rmse = line.rstrip('\n')
#        rmse_file_fptr.close()
#        rmse_sum += float(rmse)
#        rmse_sum2 += float(rmse) * float(rmse)
#
#    if nfold_exist == 0:
#        print(f"warning: nfold_exist == 0")
#    else:
#        if nfold_exist != nfold:
#            print(f"warning: nfold_exist = {nfold_exist}, nfold = {nfold}")
#        rmse_ave = rmse_sum / nfold_exist
#        rmse_var = (rmse_sum2 - rmse_sum * rmse_sum / nfold_exist) / nfold_exist
#        rmse_stddev = math.sqrt(rmse_var)
#        print(f"{mu:.1e} {rmse_ave} {rmse_stddev}")
#        print(f"{mu:.1e} {rmse_ave} {rmse_stddev}", file=mu_rmse_file_fptr)
#
#
#print("log x on", file=mu_rmse_file_fptr)
#print("lw 7", file=mu_rmse_file_fptr)
#print("line on", file=mu_rmse_file_fptr)
#print("ma 6 on", file=mu_rmse_file_fptr)
#
#mu_rmse_file_fptr.close()
#
#
##
##if __name__ == "__main__":
##    main()
#y
