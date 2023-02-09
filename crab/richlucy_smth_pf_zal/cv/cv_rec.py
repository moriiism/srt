#!/usr/bin/env python3
#
# cv_rec.py
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
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import time

def run_cmd(cmd):
    subprocess.call(cmd)

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
nem            = int(sys.argv[iarg]); iarg += 1
tol_em         = float(sys.argv[iarg]); iarg += 1
mu_list        = sys.argv[iarg]; iarg += 1
gamma_list     = sys.argv[iarg]; iarg += 1
acc_method     = sys.argv[iarg]; iarg += 1
cpu_num        = int(sys.argv[iarg]); iarg += 1

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
print("nem = ", nem)
print("tol_em = ", tol_em)
print("mu_list = ", mu_list)
print("gamma_list = ", gamma_list)
print("acc_method = ", acc_method)
print("cpu_num = ", cpu_num)

cpu_num_recommend_max = os.cpu_count() // 2
if(cpu_num > cpu_num_recommend_max):
    exit()
    
srt_dir = "/home/morii/work/github/moriiism/srt"

# data_list
print("data_list ...")
data_file_lst = []
phase_id_lst = []
phase_tag_lst = []
phase_ratio_lst = []
live_time_ratio_lst = []
flux_0_lst = []
flux_0_err_lst = []
data_list_fptr = open(data_list, "r")
for line in data_list_fptr:
    line = line.rstrip()
    if(line[0] == "#"):
        continue
    print(line)
    (data_file, phase_id, phase_tag,
     phase_ratio, live_time_ratio,
     flux_0, flux_0_err) = line.split()
    data_file_lst.append(data_file)
    phase_id_lst.append(phase_id)
    phase_tag_lst.append(phase_tag)
    phase_ratio_lst.append(phase_ratio)
    live_time_ratio_lst.append(live_time_ratio)
    flux_0_lst.append(flux_0)
    flux_0_err_lst.append(flux_0_err)
    
data_list_fptr.close()
print(data_file_lst)
print(phase_id_lst)
print(phase_tag_lst)
print(phase_ratio_lst)
print(live_time_ratio_lst)
print(flux_0_lst)
print(flux_0_err_lst)
print("data_list ... done.")

# make observation image for cross-validation
print("obs_cv ...")
for index in range(len(data_file_lst)):
    outdir_mkobs_cv = (outdir + "/" + "obs_cv" + "/"
                       + phase_id_lst[index] + "_"
                       + phase_tag_lst[index])
    outfile_head_mkobs_cv = "sim"
    cmd = [srt_dir + "/" + "mkobs_cv/mkobs_cv",
           data_file_lst[index], str(rand_seed), str(nfold),
           outdir_mkobs_cv, outfile_head_mkobs_cv]
    print(cmd)
    subprocess.call(cmd)
print("obs_cv ... done.")

# make response matrix, normed response matrix, 
# and efficiency matrix files
outdir_resp = outdir + "/" + "resp"
outfile_head_resp = "cv"
nphoton_input_resp = 100
cmd = [srt_dir + "/" + "mkresp/mkresp",
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
cmd = [srt_dir + "/" + "mkimg_points/mkimg_points",
       point_src_dat_file, outdir_point_src,
       outfile_head_point_src,
       str(nskyx), str(nskyy)]
print(cmd)
subprocess.call(cmd)

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


cmd = ["mkdir", outdir + "/" + "rec"]
print(cmd)
subprocess.call(cmd)

cmd_lst = []
for mu in mu_par_lst:
    print("cv_rec.py: mu = ", mu)
    mu_dir = f"mu{mu:.1e}"
    cmd = ["mkdir", outdir + "/" + "rec" + "/"
           + mu_dir]
    print(cmd)
    subprocess.call(cmd)
    for gamma in gamma_par_lst:
        print("cv_rec.py: mu = ", mu, " gamma = ", gamma)
        gamma_dir = f"gamma{gamma:.1e}"
        cmd = ["mkdir", outdir + "/" + "rec" + "/"
               + mu_dir + "/" + gamma_dir]
        print(cmd)
        subprocess.call(cmd)
        nfold_dir = f"nfold{nfold:02}"
        cmd = ["mkdir", outdir + "/" + "rec" + "/"
               + mu_dir + "/" + gamma_dir + "/"
               + nfold_dir]
        print(cmd)
        subprocess.call(cmd)
        for ifold in range(nfold):
            print("cv_rec.py: mu = ", mu, " gamma = ", gamma,
                  f" nfold({nfold:02}):ifold({ifold:02})")
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
            print(f"# data_file  data_vl_file  phase_tag  phase_ratio  live_time_ratio  flux_0  flux_0_err",
                  file=data_list_fptr)
            for index in range(len(phase_tag_lst)):
                data_file = (
                    outdir + "/" + "obs_cv" + "/"
                    + phase_id_lst[index] + "_"
                    + phase_tag_lst[index] + "/"
                    + f"sim_obs_{rand_seed:04}_"
                    + f"{nfold:02}fold{ifold:02}_tr.fits")
                data_vl_file = (
                    outdir + "/" + "obs_cv" + "/"
                    + phase_id_lst[index] + "_"
                    + phase_tag_lst[index] + "/"
                    + f"sim_obs_{rand_seed:04}_"
                    + f"{nfold:02}fold{ifold:02}_vl.fits")
                print(f"{data_file} {data_vl_file} " +
                      f"{phase_tag_lst[index]} " +
                      f"{phase_ratio_lst[index]} " +
                      f"{live_time_ratio_lst[index]} " +
                      f"{flux_0_lst[index]} " +
                      f"{flux_0_err_lst[index]}",
                      file=data_list_fptr)
            data_list_fptr.close()

            fixed_src_norm_file = (
                f"{outdir}/skyorg" + "/" +
                f"{outfile_head_point_src}_norm.fits")
            resp_norm_file = (
                f"{outdir}/resp" + "/" +
                "cv_resp_norm.fits")
            eff_file = (
                f"{outdir}/resp" + "/" +
                "cv_eff.fits")
            outdir_rec = (
                f"{outdir}/rec/mu{mu:.1e}/gamma{gamma:.1e}/" +
                f"nfold{nfold:02}/ifold{ifold:02}")
            outfile_head_rec = "rl"
            cmd = [srt_dir + "/" + "crab" + "/"
                   + "richlucy_smth_pf_zal" + "/"
                   + "richlucy_smth_pf_zal_cuda",
                   data_list, bg_file,
                   fixed_src_norm_file, resp_norm_file, eff_file,
                   str(nskyx), str(nskyy), str(ndetx), str(ndety),
                   outdir_rec, outfile_head_rec,
                   str(nem), str(tol_em),
                   str(mu), str(gamma), acc_method, str(nfold)]
            print(cmd)
            cmd_lst.append(cmd)

print("len of cmd_lst = ", len(cmd_lst))
with ProcessPoolExecutor(max_workers=cpu_num) as executor:
    #result = list(executor.map(run_cmd, cmd_lst),
    #              total=len(cmd_lst))
    executor.map(run_cmd, cmd_lst)


# reconstruct for full data
for mu in mu_par_lst:
    print("mu = ", mu)
    mu_dir = f"mu{mu:.1e}"
    for gamma in gamma_par_lst:
        print("  gamma = ", gamma)
        gamma_dir = f"gamma{gamma:.1e}"

        # make data.list
        data_list = (outdir + "/" + "rec" + "/"
                     + mu_dir + "/" + gamma_dir + "/"
                     + "data.list")
        data_list_fptr = open(data_list, "w")
        print(f"# data_file  data_vl_file  phase_tag  phase_ratio  live_time_ratio  flux0  flux0_err",
              file=data_list_fptr)
        for index in range(len(phase_tag_lst)):
            data_vl_file = data_file_lst[index];
            print(f"{data_file_lst[index]} {data_vl_file} " +
                  f"{phase_tag_lst[index]} " +
                  f"{phase_ratio_lst[index]} " +
                  f"{live_time_ratio_lst[index]} " +
                  f"{flux_0_lst[index]} " +
                  f"{flux_0_err_lst[index]}",
                  file=data_list_fptr)
        data_list_fptr.close()
        fixed_src_norm_file = (
            f"{outdir}/skyorg" + "/" +
            f"{outfile_head_point_src}_norm.fits")
        resp_norm_file = (
            f"{outdir}/resp" + "/" +
            "cv_resp_norm.fits")
        eff_file = (
            f"{outdir}/resp" + "/" +
            "cv_eff.fits")
        outdir_rec = (
            f"{outdir}/rec/mu{mu:.1e}/gamma{gamma:.1e}")
        outfile_head_rec = "rl"
        nfold = 0
        cmd = [srt_dir + "/" + "crab" + "/"
               + "richlucy_smth_pf_zal" + "/"
               + "richlucy_smth_pf_zal_cuda",
               data_list, bg_file,
               fixed_src_norm_file, resp_norm_file, eff_file,
               str(nskyx), str(nskyy), str(ndetx), str(ndety),
               outdir_rec, outfile_head_rec,
               str(nem), str(tol_em),
               str(mu), str(gamma), acc_method, str(nfold)]
        print(cmd)
        subprocess.call(cmd)
        

#if __name__ == "__main__":
#    main()

