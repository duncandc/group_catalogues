#!/bin/bash

python make_valueadded_yang_groupcat_2.py sample1_L_model &
python make_valueadded_yang_groupcat_2.py sample2_L_model &
python make_valueadded_yang_groupcat_2.py sample3_L_model &

python make_valueadded_yang_groupcat_1.py sample1_M_model &
python make_valueadded_yang_groupcat_1.py sample2_M_model &
python make_valueadded_yang_groupcat_1.py sample3_M_model &