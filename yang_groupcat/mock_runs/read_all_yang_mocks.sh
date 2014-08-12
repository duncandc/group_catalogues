#!/bin/bash

#version 1
python read_yang_mockruns_groupcat.py 0 1 &
python read_yang_mockruns_groupcat.py 1 1 &

#version 2
python read_yang_mockruns_groupcat.py 0 2 &
python read_yang_mockruns_groupcat.py 1 2 &

#version 3
python read_yang_mockruns_groupcat.py 0 3 &
python read_yang_mockruns_groupcat.py 1 3 &

#version 4
python read_yang_mockruns_groupcat.py 0 4 &
python read_yang_mockruns_groupcat.py 1 4 &

#version 5
python read_yang_mockruns_groupcat.py 0 5 &
python read_yang_mockruns_groupcat.py 1 5 &

