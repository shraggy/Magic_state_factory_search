#!/bin/zsh
cd "/Users/shraddhasingh/Library/CloudStorage/OneDrive-IBM/Documents/GitHub_shraggy/Magic_state_factory_search"
nice -n 19 python3 searches/staged.py --L 3 --k 1 2 3 4 --out outputs/staged --jobs 3 > outputs/staged/L0_l3.log 2>&1 &
nice -n 19 python3 searches/staged.py --L 2 --k 1 2 3 4 5 6 --out outputs/staged --jobs 3 > outputs/staged/L0_l2.log 2>&1 &
nice -n 19 python3 searches/staged.py --L 4 --k 1 2 3 --out outputs/staged --jobs 3 > outputs/staged/L0_l4.log 2>&1 &
wait
echo L0 done
