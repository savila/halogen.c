reset


set term postscript enhance
set output "MD_alpha.eps"

set xrange [0:200]
set multiplot layout 2,1

#set xrange [0:200]
#set xrange [15:120]
set yrange [-20:80]

set ylabel 'r^2 {/Symbol x} (r)'

#set title '2pt CF in real space'

set tmargin at screen  0.9
set bmargin at screen  0.5
set lmargin at screen  0.1
set rmargin at screen  0.95



plot "/bigdata/savila/BIG_CATS/MultiDark/MultiDark_1Gp/FoF_n3.5e-4_halos_lin_2PCF" u 1:($2*$1*$1) w points ps 1.5 lc 1 lw 2 t 'MD FoF' , \
"/bigdata/savila/BIG_CATS/MultiDark/FoF_BigMD_n3.5e-4.albert_2PCF" u 1:($2*$1*$1) w points ps 1 lc 2 lw 2 t 'BigMD FoF' , \
"output/MD_FOF_l4_alphaBigMD_run0.2PCF" u 1:($2*$1*$1) w lines lt 1 lc 2 t 'sim-run1 alpha-BMD', \
"output/MD_FOF_l4_N512_run0.2PCF" u 1:($2*$1*$1) w lines lt 1 lc 1 t 'sim-N512 alpha-N512', \
"output/MD_FOF_l4_mean.2PCF" u 1:($2*$1*$1) w lines lt 2 lc 3 t 'sim-run1 alpha-run1', \
"output/MDFOFrun2_l4_run0.2PCF" u 1:($2*$1*$1) w lines lt 2 lc 4 t 'sim-run2 alpha-run1', \
"output/MDFOFrun3_l4_run0.2PCF" u 1:($2*$1*$1) w lines lt 2 lc 5 t 'sim-run3 alpha-run1', \
"/home/savila/HALOGEN/c_halogen/v0.4/output/single_alpha1.30_L1000_N512_0.74_NC64.2PCF" u 1:($2*$1*$1) w lines lt 3 lc 6 t 'v0.4', \
"/home/savila/HALOGEN/c_halogen/v0.4/output/single_alpha1.29_L1000_N512_0.74_NC64.2PCF" u 1:($2*$1*$1) w lines lt 3 lc 7 t 'v0.4 1.29'


#"output/MD_FOF_l7_mean.2PCF" u 1:($2*$1*$1) w lines lt 1 lc 5 t 'l=7', \
#"output/check_v0.7.5_l3_run0.2PCF" u 1:($2*$1*$1) w lines lt 3 lc 10 t 'l=3 old'

set ylabel 'ratio'
set xlabel 'r'

set yrange [0.8:1.2]
#set yrange [0.0:2.0]


set tmargin at screen  0.45
set bmargin at screen  0.1
set lmargin at screen  0.1
set rmargin at screen  0.95


plot "<paste /bigdata/savila/BIG_CATS/MultiDark/FoF_BigMD_n3.5e-4.albert_2PCF output/MD_FOF_l4_alphaBigMD_run0.2PCF" u 1:($6/$2) w lines lc 2 lt 1 lw 1  t '', 1.0 w lines lt 1 lc -1 t '', \
"<paste /bigdata/savila/BIG_CATS/MultiDark/MultiDark_1Gp/FoF_n3.5e-4_halos_lin_2PCF output/MD_FOF_l4_N512_run0.2PCF" u 1:($6/$2) w lines lc 1 lt 1 lw 1 t '', \
"<paste /bigdata/savila/BIG_CATS/MultiDark/MultiDark_1Gp/FoF_n3.5e-4_halos_lin_2PCF output/MD_FOF_l4_mean.2PCF" u 1:($6/$2) w lines lc 3 lt 2 lw 1 t '', \
"<paste /bigdata/savila/BIG_CATS/MultiDark/MultiDark_1Gp/FoF_n3.5e-4_halos_lin_2PCF output/MDFOFrun2_l4_run0.2PCF" u 1:($6/$2) w lines lc 4 lt 2 lw 1 t '', \
"<paste /bigdata/savila/BIG_CATS/MultiDark/MultiDark_1Gp/FoF_n3.5e-4_halos_lin_2PCF output/MDFOFrun3_l4_run0.2PCF" u 1:($6/$2) w lines lc 5 lt 2 lw 1 t ''

#plot "<paste /bigdata/savila/BIG_CATS/MultiDark/MultiDark_1Gp/FoF_n3.5e-4_halos_lin_2PCF output/MD_FOF_l4_alphaBigMD_run0.2PCF" u 1:($6/$2) w lines lc 2 lt 1 lw 1  t '', 1.0 w lines lt 1 lc -1 t '', \

