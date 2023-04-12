#!/usr/bin/env python3
#
# run_cv.py
#
'''
Overview:
  run cv_rec.py and cv_smr.py for many data
Input:

Output:

Data:

Details:

'''
import os
import sys
import subprocess
import math
import numpy as np

srt_dir = "/home/morii/work/github/moriiism/srt"

def LoadList1(list_file):
    line_lst = []
    list_fptr = open(list_file, "r")
    for line in list_fptr:
        line = line.rstrip()
        if(line == ""):
            continue
        if(line[0] == "#"):
            continue
        # print(line)
        line_lst.append(line)

    list_fptr.close()
    return(line_lst)

def LoadList2(list_file):
    line1_lst = []
    line2_lst = []
    list_fptr = open(list_file, "r")
    for line in list_fptr:
        line = line.rstrip()
        if(line == ""):
            continue
        if(line[0] == "#"):
            continue
        (line1, line2) = line.split()
        line1_lst.append(line1)
        line2_lst.append(line2)

    list_fptr.close()
    return(line1_lst, line2_lst)

def LoadLiveTimeRatioList(live_time_ratio_list, phase_id_lst):
    (phase_id_ltr_lst,
     live_time_ratio_lst) = LoadList2(live_time_ratio_list)
    live_time_ratio_out_lst = []
    for phase_id in phase_id_lst:
        for index in range(len(phase_id_ltr_lst)):
            if (phase_id == phase_id_ltr_lst[index]):
                live_time_ratio_out_lst.append(
                    live_time_ratio_lst[index])
    return(live_time_ratio_out_lst)


def LoadPhaseListFile(phase_list_file):
    phase_id_lst = []
    phase_st_lst = []
    phase_ed_lst = []
    phase_ratio_lst = []
    phase_tag_lst = []    
    phase_list_fptr = open(phase_list_file, "r")
    for line in phase_list_fptr:
        line = line.rstrip()
        if(line == ""):
            continue
        if(line[0] == "#"):
            continue
        print(line)
        (phase_id, phase_st, phase_ed,
         phase_ratio, phase_tag) = line.split()
        phase_id_lst.append(phase_id)
        phase_st_lst.append(phase_st)
        phase_ed_lst.append(phase_ed)
        phase_ratio_lst.append(phase_ratio)        
        phase_tag_lst.append(phase_tag)
    
    phase_list_fptr.close()
    print(phase_id_lst)
    print(phase_st_lst)
    print(phase_ed_lst)
    print(phase_ratio_lst)    
    print(phase_tag_lst)

    return (phase_id_lst, phase_st_lst, phase_ed_lst,
            phase_ratio_lst, phase_tag_lst)


#
# main
#

iarg = 1
work_dir       = sys.argv[iarg]; iarg += 1
target_list    = sys.argv[iarg]; iarg += 1
data_ver_list  = sys.argv[iarg]; iarg += 1
resp_ver_list  = sys.argv[iarg]; iarg += 1
ene_band_list  = sys.argv[iarg]; iarg += 1
rec_prog_list  = sys.argv[iarg]; iarg += 1
mu_list        = sys.argv[iarg]; iarg += 1
gamma_list     = sys.argv[iarg]; iarg += 1
phase_list     = sys.argv[iarg]; iarg += 1
nfold          = int(sys.argv[iarg]); iarg += 1
nskyx          = int(sys.argv[iarg]); iarg += 1
nskyy          = int(sys.argv[iarg]); iarg += 1
ndetx          = int(sys.argv[iarg]); iarg += 1
ndety          = int(sys.argv[iarg]); iarg += 1
use_cuda       = int(sys.argv[iarg]); iarg += 1

print("work_dir = ", work_dir)
print("target_list = ", target_list)
print("data_ver_list = ", data_ver_list)
print("resp_ver_list = ", resp_ver_list)
print("ene_band_list = ", ene_band_list)
print("rec_prog_list = ", rec_prog_list)
print("mu_list = ", mu_list)
print("gamma_list = ", gamma_list)
print("phase_list = ", phase_list)
print("nfold = ", nfold)
print("nskyx = ", nskyx)
print("nskyy = ", nskyy)
print("ndetx = ", ndetx)
print("ndety = ", ndety)
print("use_cuda = ", use_cuda)

# target object
# data version
# resp version
# energy band
# reconstruction program
target_lst = LoadList1(target_list)
data_ver_lst = LoadList1(data_ver_list)
resp_ver_lst = LoadList1(resp_ver_list)
ene_band_lst = LoadList1(ene_band_list)
rec_prog_lst = LoadList1(rec_prog_list)
mu_lst = LoadList1(mu_list)
gamma_lst = LoadList1(gamma_list)

(phase_id_lst, phase_st_lst, phase_ed_lst,
 phase_ratio_lst, phase_tag_lst) = LoadPhaseListFile(phase_list)

outdir_lst = []
target_joint_lst = []
data_ver_joint_lst = []
resp_ver_joint_lst = []
ene_band_joint_lst = []
rec_prog_joint_lst = []
for target in target_lst:
        for data_ver in data_ver_lst:
            for resp_ver in resp_ver_lst:        
                for ene_band in ene_band_lst:
                    for rec_prog in rec_prog_lst:
                        outdir = (work_dir + "/"
                                   + target + "/"
                                   + "hxt12" + "/"
                                   + "data_" + data_ver + "/"
                                   + "resp_" + resp_ver + "/"
                                   + ene_band + "/"
                                   + rec_prog)
                        outdir_lst.append(outdir)
                        target_joint_lst.append(target)
                        data_ver_joint_lst.append(data_ver)
                        resp_ver_joint_lst.append(resp_ver)
                        ene_band_joint_lst.append(ene_band)
                        rec_prog_joint_lst.append(rec_prog)

for outdir in outdir_lst:
    cmd = ["mkdir", "-p", outdir]
    print(cmd)
    subprocess.call(cmd)

# ymaeda share data
ymaeda_share_dir = "/home/ymaeda/work/share/crab_hxt"

for index in range(len(outdir_lst)):
    print(outdir_lst[index])
    data_ver = data_ver_joint_lst[index]
    resp_ver = resp_ver_joint_lst[index]    
    ene_band = ene_band_joint_lst[index]
    live_time_ratio_det1_list = (ymaeda_share_dir + "/"
                                 + "data" + "/"
                                 + data_ver + "/"
                                 + "hxt1" + "/"
                                 + "live_time_ratio.list")
    live_time_ratio_det1_lst = LoadLiveTimeRatioList(
        live_time_ratio_det1_list,
        phase_id_lst)

    live_time_ratio_det2_list = (ymaeda_share_dir + "/"
                                 + "data" + "/"
                                 + data_ver + "/"
                                 + "hxt2" + "/"
                                 + "live_time_ratio.list")
    live_time_ratio_det2_lst = LoadLiveTimeRatioList(
        live_time_ratio_det2_list,
        phase_id_lst)
    
    # find phase_off
    flux_det1_lst = []
    flux_det2_lst = []
    flux_err_det1_lst = []
    flux_err_det2_lst = []
    for index_phase in range(len(phase_id_lst)):
        phase_ratio = phase_ratio_lst[index_phase]
        live_time_ratio_det1 = live_time_ratio_det1_lst[index_phase]
        # hxt1_phase_42_000_2047_x1178y1175_w80.fits
        data_on_file = (ymaeda_share_dir + "/"
                        + "data" + "/"
                        + data_ver + "/"
                        + "hxt1" + "/"
                        + "hxt1" + "_"
                        + "phase" + "_"
                        + phase_id_lst[index_phase] + "_"
                        + ene_band + "_"
                        + "x1178y1175_w80.fits")
        cmd = [srt_dir + "/" + "crab" + "/"
               + "calc_flux0" + "/"
               + "calc_flux0",
               data_on_file,
               phase_ratio,
               live_time_ratio_det1,
               "none", str(0.0), str(0.0)]
        print(cmd)
        # subprocess.call(cmd)
        output_str = subprocess.run(
            cmd, capture_output=True, text=True).stdout
        
        (flux_det1, flux_err_det1) = (
            output_str.splitlines()[-1].split())
        flux_det1_lst.append(flux_det1)
        flux_err_det1_lst.append(flux_err_det1)

        # hxt2
        live_time_ratio_det2 = live_time_ratio_det2_lst[index_phase]
        # hxt2_phase_42_000_2047_x1178y1175_w80.fits
        data_on_file = (ymaeda_share_dir + "/"
                        + "data" + "/"
                        + data_ver + "/"
                        + "hxt2" + "/"
                        + "hxt2" + "_"
                        + "phase" + "_"
                        + phase_id_lst[index_phase] + "_"
                        + ene_band + "_"
                        + "x1178y1175_w80.fits")
        cmd = [srt_dir + "/" + "crab" + "/"
               + "calc_flux0" + "/"
               + "calc_flux0",
               data_on_file,
               phase_ratio,
               live_time_ratio_det2,
               "none", str(0.0), str(0.0)]
        print(cmd)
        # subprocess.call(cmd)
        output_str = subprocess.run(
            cmd, capture_output=True, text=True).stdout
        (flux_det2, flux_err_det2) = (
            output_str.splitlines()[-1].split())
        flux_det2_lst.append(flux_det2)
        flux_err_det2_lst.append(flux_err_det2)

        
    # (flux_det1 + flux_det2) / 2
    flux_lst = []
    for ind in range(len(flux_det1_lst)):
        flux_ave = (float(flux_det1_lst[ind]) +
                    float(flux_det2_lst[ind])) / 2.0
        flux_lst.append(flux_ave)

    print(flux_lst)
    index_off = flux_lst.index(min(flux_lst))
    index_max = flux_lst.index(max(flux_lst))
    print("index_off = ", index_off)
    print("index_max = ", index_max)

    # calc flux0_lst, flux0_err_lst:
    # flux0 = flux_on - flux_off
    flux0_det1_lst = []
    flux0_det2_lst = []    
    flux0_err_det1_lst = []
    flux0_err_det2_lst = []    
    for index_phase in range(len(phase_id_lst)):
        phase_ratio_on = phase_ratio_lst[index_phase]
        live_time_ratio_on_det1 = (
            live_time_ratio_det1_lst[index_phase])
        phase_ratio_off = phase_ratio_lst[index_off]
        live_time_ratio_off_det1 = (
            live_time_ratio_det1_lst[index_off])
        data_on_file = (ymaeda_share_dir + "/"
                        + "data" + "/"
                        + data_ver + "/"
                        + "hxt1" + "/"
                        + "hxt1" + "_"
                        + "phase" + "_"
                        + phase_id_lst[index_phase] + "_"
                        + ene_band + "_"
                        + "x1178y1175_w80.fits")
        data_off_file = (ymaeda_share_dir + "/"
                         + "data" + "/"
                         + data_ver + "/"
                         + "hxt1" + "/"
                         + "hxt1" + "_"
                         + "phase" + "_"
                         + phase_id_lst[index_off] + "_"
                         + ene_band + "_"
                         + "x1178y1175_w80.fits")
        cmd = [srt_dir + "/" + "crab" + "/"
               + "calc_flux0" + "/"
               + "calc_flux0",
               data_on_file,
               phase_ratio_on,
               live_time_ratio_on_det1,
               data_off_file,
               phase_ratio_off,
               live_time_ratio_off_det1]
        print(cmd)
        # subprocess.call(cmd)
        output_str = subprocess.run(
            cmd, capture_output=True, text=True).stdout
        (flux0_det1, flux0_err_det1) = (
            output_str.splitlines()[-1].split())
        flux0_det1_lst.append(flux0_det1)
        flux0_err_det1_lst.append(flux0_err_det1)

        # hxt2
        phase_ratio_on = phase_ratio_lst[index_phase]
        live_time_ratio_on_det2 = (
            live_time_ratio_det2_lst[index_phase])
        phase_ratio_off = phase_ratio_lst[index_off]
        live_time_ratio_off_det2 = (
            live_time_ratio_det2_lst[index_off])
        data_on_file = (ymaeda_share_dir + "/"
                        + "data" + "/"
                        + data_ver + "/"
                        + "hxt2" + "/"
                        + "hxt2" + "_"
                        + "phase" + "_"
                        + phase_id_lst[index_phase] + "_"
                        + ene_band + "_"
                        + "x1178y1175_w80.fits")
        data_off_file = (ymaeda_share_dir + "/"
                         + "data" + "/"
                         + data_ver + "/"
                         + "hxt2" + "/"
                         + "hxt2" + "_"
                         + "phase" + "_"
                         + phase_id_lst[index_off] + "_"
                         + ene_band + "_"
                         + "x1178y1175_w80.fits")
        cmd = [srt_dir + "/" + "crab" + "/"
               + "calc_flux0" + "/"
               + "calc_flux0",
               data_on_file,
               phase_ratio_on,
               live_time_ratio_on_det2,
               data_off_file,
               phase_ratio_off,
               live_time_ratio_off_det2]
        print(cmd)
        # subprocess.call(cmd)
        output_str = subprocess.run(
            cmd, capture_output=True, text=True).stdout
        (flux0_det2, flux0_err_det2) = (
            output_str.splitlines()[-1].split())
        flux0_det2_lst.append(flux0_det2)
        flux0_err_det2_lst.append(flux0_err_det2)

    print(flux0_det1_lst)
    print(flux0_err_det1_lst)    
    print(flux0_det2_lst)
    print(flux0_err_det2_lst)

    # (flux0_det1 + flux0_det2) / 2
    flux0_lst = []
    flux0_err_lst = []
    for ind in range(len(flux0_det1_lst)):
        flux0_ave = (float(flux0_det1_lst[ind]) +
                     float(flux0_det2_lst[ind])) / 2.0
        flux0_ave_err = sqrt(
            float(flux0_err_det1_lst[ind]) *
            float(flux0_err_det1_lst[ind]) +
            float(flux0_err_det2_lst[ind]) *
            float(flux0_err_det2_lst[ind])) / 2.0
        
        flux0_lst.append(flux0_ave)
        flux0_err_lst.append(flux0_ave_err)
    
    #
    #-----> ds9: (51, 51) --> (50, 50)
    #
    print(outdir_lst[index])
    
    data_list = (outdir_lst[index] + "/" + "data.list")
    data_list_fptr = open(data_list, "w")
    print("# data_file1  data_file2  phase_id  phase_tag  " +
          "phase_ratio  " +
          "live_time_ratio_det1  live_time_ratio_det2  " +
          "flux0  flux0_err",
          file=data_list_fptr)
    for index_phase in range(len(phase_id_lst)):
        data_file1 = (ymaeda_share_dir + "/"
                      + "data" + "/"
                      + data_ver + "/"
                      + "hxt1" + "/"
                      + "hxt1" + "_"
                      + "phase" + "_"
                      + phase_id_lst[index_phase] + "_"
                      + ene_band + "_"
                      + "x1178y1175_w80.fits")
        data_file2 = (ymaeda_share_dir + "/"
                      + "data" + "/"
                      + data_ver + "/"
                      + "hxt2" + "/"
                      + "hxt2" + "_"
                      + "phase" + "_"
                      + phase_id_lst[index_phase] + "_"
                      + ene_band + "_"
                      + "x1178y1175_w80.fits")
        print(f"{data_file1}  {data_file2}  " +
              f"{phase_id_lst[index_phase]} " +
              f"{phase_tag_lst[index_phase]} " +
              f"{phase_ratio_lst[index_phase]} " +
              f"{live_time_ratio_det1_lst[index_phase]} " +
              f"{live_time_ratio_det2_lst[index_phase]} " +
              f"{flux0_lst[index_phase]} " +
              f"{flux0_err_lst[index_phase]}",
              file=data_list_fptr)
    data_list_fptr.close()

    bg_file1 = "none"
    bg_file2 = "none"
    rand_seed = 1
    respdir1 = (ymaeda_share_dir + "/" + "resp" + "/"
                + resp_ver + "/" + "hxt1" + "/"
                + ene_band + "/" + "gimage_100")
    respdir2 = (ymaeda_share_dir + "/" + "resp" + "/"
                + resp_ver + "/" + "hxt2" + "/"
                + ene_band + "/" + "gimage_100")
    posx_point_src = 50
    posy_point_src = 50
    outdir_rl = (outdir_lst[index] + "/" + "cv")
    # nem = 10000
    nem = 1000
    tol_em = 1.0e-7
    acc_method = "none"
    cpu_num = 1

    cmd = [srt_dir + "/" + "crab" + "/" +
           "richlucy_smth_pf_zal_det2" + "/" +
           "cv" + "/" + "cv_rec.py",
           data_list, bg_file1, bg_file2,
           str(rand_seed), str(nfold),
           str(nskyx), str(nskyy), str(ndetx), str(ndety),
           respdir1, respdir2,
           str(posx_point_src), str(posy_point_src),
           outdir_rl, str(nem), str(tol_em),
           mu_list, gamma_list,
           acc_method, str(cpu_num), str(use_cuda)]
    print(cmd)
    subprocess.call(cmd)

    cmd = [srt_dir + "/" + "crab" + "/" +
           "richlucy_smth_pf_zal_det2" + "/" +
           "cv" + "/" + "cv_smr.py",
           str(nfold),
           outdir_rl,
           mu_list, gamma_list]
    print(cmd)
    subprocess.call(cmd)


    
#
#if __name__ == "__main__":
#    main()

