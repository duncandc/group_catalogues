#!/bin/bash

python make_valueadded_yang_groupcat_L.py sample1_L_model &
python make_valueadded_yang_groupcat_L.py sample2_L_model &
python make_valueadded_yang_groupcat_L.py sample3_L_model &

python make_valueadded_yang_groupcat_M.py sample1_M_model &
python make_valueadded_yang_groupcat_M.py sample2_M_model &
python make_valueadded_yang_groupcat_M.py sample3_M_model &
