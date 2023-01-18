#!/usr/bin/env python3
#
# cv.py
#
'''
Overview:

Input:

Output:

Data:

Details:

Example:
 

'''
import os
import sys
import subprocess
#import pandas as pd
#import math

iarg = 1
data_list      = sys.argv[iarg]; iarg += 1
bg_file        = sys.argv[iarg]; iarg += 1
rand_seed      = int(sys.argv[iarg]); iarg += 1
nfold          = int(sys.argv[iarg]); iarg += 1
nskyx          = int(sys.argv[iarg]); iarg += 1
nskyy          = int(sys.argv[iarg]); iarg += 1
ndetx          = int(sys.argv[iarg]); iarg += 1
ndety          = int(sys.argv[iarg]); iarg += 1
respdir        = sys.argv[iarg]; iarg += 1
posx_point_src = int(sys.argv[iarg]); iarg += 1
posy_point_src = int(sys.argv[iarg]); iarg += 1
outdir         = sys.argv[iarg]; iarg += 1
outfile_head   = sys.argv[iarg]; iarg += 1
nem            = int(sys.argv[iarg]); iarg += 1
tol_em         = float(sys.argv[iarg]); iarg += 1
mu_list        = sys.argv[iarg]; iarg += 1
gamma_list     = sys.argv[iarg]; iarg += 1
acc_method     = sys.argv[iarg]; iarg += 1

print("data_list = ", data_list)
print("bg_file = ", bg_file)
print("rand_seed = ", rand_seed)
print("nfold = ", nfold)
print("nskyx = ", nskyx)
print("nskyy = ", nskyy)
print("ndetx = ", ndetx)
print("ndety = ", ndety)
print("respdir = ", respdir)
print("posx_point_src = ", posx_point_src)
print("posy_point_src = ", posy_point_src)
print("outdir = ", outdir)
print("outfile_head = ", outfile_head)
print("nem = ", nem)
print("tol_em = ", tol_em)
print("mu_list = ", mu_list)
print("gamma_list = ", gamma_list)
print("acc_method = ", acc_method)

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

# make observation image for cross-validation
print("obs_cv ...")
for index in range(len(data_file_lst)):
    outdir_mkobs_cv = outdir + "/" + "obs_cv" + "/" + phase_tag_lst[index]
    outfile_head_mkobs_cv = outfile_head
    print(outdir_mkobs_cv)
    cmd = [srt_dir + "/" + "mkobs_cv/mkobs_cv",
           data_file_lst[index], str(rand_seed), str(nfold),
           outdir_mkobs_cv, outfile_head_mkobs_cv]
    print(cmd)
    subprocess.call(cmd)

print("obs_cv ... done.")


# make response matrix, normed response matrix, 
# and efficiency matrix files
outdir_resp = outdir + "/" + "resp"
outfile_head_resp = "hxt1"
nphoton_input_resp = 100
cmd = ["/home/morii/work/github/moriiism/srt/mkresp/mkresp", 
       respdir, outdir_resp, outfile_head_resp,
       str(nskyx), str(nskyy), str(nphoton_input_resp)]
print(cmd)
subprocess.call(cmd)


# mkimg_points
cmd = ["mkdir", outdir + "/" + "skyorg"]
print(cmd)
subprocess.call(cmd)
point_src_dat_file = (outdir + "/" + "skyorg" + "/"
                      + "crab_pulsar.dat")
point_src_dat_file_fptr = open(point_src_dat_file, "w")
print(f"{posx_point_src} {posy_point_src} 1.0",
      file=point_src_dat_file_fptr)
point_src_dat_file_fptr.close()

outdir_point_src = outdir + "/" + "skyorg"
outfile_head_point_src = "crab_pulsar"
cmd = ["/home/morii/work/github/moriiism/srt/mkimg_points/mkimg_points",
       point_src_dat_file, outdir_point_src,
       outfile_head_point_src,
       str(nskyx), str(nskyy)]
print(cmd)
subprocess.call(cmd)

# hyper_par_list
#hyper_par_lst = []
#hyper_par_fptr = open(hyper_par_list, "r")
#for line in hyper_par_fptr:
#    mu = line.rstrip('\n')
#    mu_lst.append(mu)
#mu_list_file_fptr.close()


mu_lst = [1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0]
gamma_lst = [1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1, 1.0]

#mu_lst = [1.0e-3]
#gamma_lst = [1.0e-3]

cmd = ["mkdir", outdir + "/" + "rec"]
print(cmd)
subprocess.call(cmd)

for mu in mu_lst:
    print("mu = %e", mu)
    mu_dir = f"mu{mu:.1e}"
    cmd = ["mkdir", outdir + "/" + "rec" + "/"
           + mu_dir]
    print(cmd)
    subprocess.call(cmd)
    for gamma in gamma_lst:
        print("  gamma = %s", gamma)
        gamma_dir = f"gamma{gamma:.1e}"
        cmd = ["mkdir", outdir + "/" + "rec" + "/"
               + mu_dir + "/" + gamma_dir]
        print(cmd)
        subprocess.call(cmd)
        mu = float(mu)
        gamma = float(gamma)

        nfold_dir = f"nfold{nfold:02}"
        cmd = ["mkdir", outdir + "/" + "rec" + "/"
               + mu_dir + "/" + gamma_dir + "/"
               + nfold_dir]
        print(cmd)
        subprocess.call(cmd)
        for ifold in range(nfold):
            print("      nfold(%d):ifold(%d)",
                  nfold, ifold)
            ifold_dir = f"ifold{ifold:02}"
            cmd = ["mkdir", outdir + "/" + "rec" + "/"
                   + mu_dir + "/" + gamma_dir + "/"
                   + nfold_dir + "/" + ifold_dir]
            print(cmd)
            subprocess.call(cmd)
            # make data.list
            data_list = (outdir + "/" + "rec" + "/"
                         + mu_dir + "/" + gamma_dir + "/"
                         + f"nfold{nfold:02}" + "/"
                         + f"ifold{ifold:02}" + "/"
                         + "data.list")
            data_list_fptr = open(data_list, "w")
            print(f"# data_file  phase_tag  phase_ratio  flux_0",
                  file=data_list_fptr)
            for index in range(len(phase_tag_lst)):
                data_file = (
                    outdir + "/" + "obs_cv" + "/"
                    + phase_tag_lst[index] + "/"
                    + f"{outfile_head}_obs_{rand_seed:04}_"
                    + f"{nfold:02}fold{ifold:02}_tr.fits")
                print(f"{data_file} {phase_tag_lst[index]} "
                      f"{phase_ratio_lst[index]} {flux_0_lst[index]}",
                      file=data_list_fptr)

            data_list_fptr.close()

            bg_file = "none"
            fixed_src_norm_file = (
                f"{outdir}/skyorg/{outfile_head_point_src}_norm.fits")
            resp_file = (
                f"{outdir}/resp/{outfile_head_resp}_resp.fits")
            eff_file = (
                f"{outdir}/resp/{outfile_head_resp}_eff.fits")
            outdir_rec = (
                f"{outdir}/rec/mu{mu:.1e}/gamma{gamma:.1e}/" +
                f"nfold{nfold:02}/ifold{ifold:02}")
            outfile_head_rec = (f"rec")

            cmd = [srt_dir + "/" + "crab" + "/"
                   + "richlucy_smth_pf_zal" + "/"
                   + "richlucy_smth_pf_zal_cuda",
                   data_list, bg_file,
                   fixed_src_norm_file, resp_file, eff_file,
                   str(nskyx), str(nskyy), str(ndetx), str(ndety),
                   outdir_rec, outfile_head_rec,
                   str(nem), str(tol_em),
                   str(mu), str(gamma), acc_method]
            print(cmd)
            subprocess.call(cmd)

#if __name__ == "__main__":
#    main()

