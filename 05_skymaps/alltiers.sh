#!/bin/bash

input='./anti_fits/'
output='./output/'
fit='--fit 1'
#range='-M 12 -m -12'
#range='-M 3.5 -m -3.5'
range='-M 1 -m -1'
python 1d_proj.py $input"combined_t1_iteration20.fits.gz" --relerr $input"significance_t1_iteration20.fits.gz" -o $output"1d_t1_antiflat.png" $fit $range -d -90 -D -35 -z -S 3 --title "Tier 1: 310 TeV 2011-14"
#python 1d_proj.py $input"combined_t2_iteration20.fits.gz" --relerr $input"significance_t2_iteration20.fits.gz" -o $output"1d_t2_solardi.png" $fit $range -d -90 -D -35 -z -S 3 --title "Tier 2: 1.1 PeV 2011-14"
#python 1d_proj.py $input"combined_t3_iteration08.fits.gz" --relerr $input"significance_t3_iteration08.fits.gz" -o $output"1d_t3_solardi.png" $fit $range -d -90 -D -35 -z -S 3 --title "Tier 3: 2.4 PeV 2011-21"
#python 1d_proj.py $input"combined_t4_iteration20.fits.gz" --relerr $input"significance_t4_iteration20.fits.gz" -o $output"1d_t4_solardi.png" $fit $range -d -90 -D -35 -z -S 3 --title "Tier 4: 6.6 PeV 2011-21"

