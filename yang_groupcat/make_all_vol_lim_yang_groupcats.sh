#!/bin/bash

python make_mr_vol_lim_yang_groupcat.py -20.0 &
python make_mr_vol_lim_yang_groupcat.py -19.0 &
python make_mr_vol_lim_yang_groupcat.py -18.0 &

wait

python make_sm_vol_lim_yang_groupcat.py 10.6 -20.0 &
python make_sm_vol_lim_yang_groupcat.py 10.1 -19.0 &
python make_sm_vol_lim_yang_groupcat.py 9.7 -18.0 &