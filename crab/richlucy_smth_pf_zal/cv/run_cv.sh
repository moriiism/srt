#!/bin/sh

export LANG=C
source ~/work/github/moriiism/mitool/setup/setup_arb01.sh

work_dir=/home/morii/work/arb/ana/run
cd $work_dir
mkdir conf

cat << EOF > conf/target.list
crab
EOF

cat << EOF > conf/telescope.list
hxt1
hxt2
EOF

cat << EOF > conf/data_ver.list
nominal
EOF

cat << EOF > conf/resp_ver.list
nominal
m0p1
p0p1
EOF

cat << EOF > conf/ene_band.list
#000_2047
036_150
#036_300
151_300
301_700
#649_1150
EOF

cat << EOF > conf/rec_prog.list
crab_rl_smth_pf_zal
EOF

cat << EOF > conf/mu.list
1.0000e-10
1.0000e-09
1.0000e-08
1.0000e-07
1.0000e-06
1.0000e-05
1.0000e-04
1.0000e-03
1.0000e-02
1.0000e-01
1.0000e+00
EOF

cat << EOF > conf/gamma.list
1.0000e-10
1.0000e-09
1.0000e-08
1.0000e-07
1.0000e-06
1.0000e-05
1.0000e-04
1.0000e-03
1.0000e-02
1.0000e-01
1.0000e+00
EOF

# phase_list
cat << EOF > conf/phase_20.list
# phase_id phase_st phase_ed phase_ratio phase_tag
21  0.854  1.104  0.25  on2
22  0.104  0.504  0.40  off1
23  0.504  0.704  0.20  on1
24  0.704  0.854  0.15  off2
EOF

cat << EOF > conf/phase_40.list
# phase_id phase_st phase_ed phase_ratio phase_tag
40  0.000  0.100  0.10  0.0-0.1
41  0.100  0.200  0.10  0.1-0.2
42  0.200  0.300  0.10  0.2-0.3
43  0.300  0.400  0.10  0.3-0.4
44  0.400  0.500  0.10  0.4-0.5
45  0.500  0.600  0.10  0.5-0.6
46  0.600  0.700  0.10  0.6-0.7
47  0.700  0.800  0.10  0.7-0.8
48  0.800  0.900  0.10  0.8-0.9
49  0.900  1.000  0.10  0.9-1.0
EOF

work_dir=$work_dir
target_list=conf/target.list
telescope_list=conf/telescope.list
data_ver_list=conf/data_ver.list
resp_ver_list=conf/resp_ver.list
ene_band_list=conf/ene_band.list
rec_prog_list=conf/rec_prog.list
mu_list=conf/mu.list
gamma_list=conf/gamma.list
phase_list=conf/phase_20.list
#phase_list=conf/phase_40.list
nfold=5
nskyx=101
nskyy=101
ndetx=80
ndety=80

srt_dir=/home/morii/work/github/moriiism/srt
$srt_dir/crab/richlucy_smth_pf_zal/cv/run_cv.py \
$work_dir \
$target_list \
$telescope_list \
$data_ver_list \
$resp_ver_list \
$ene_band_list \
$rec_prog_list \
$mu_list \
$gamma_list \
$phase_list \
$nfold \
$nskyx \
$nskyy \
$ndetx \
$ndety
