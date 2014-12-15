what do all these scripts do?

read_berlind_groupcat_mockruns_x.x.py
    reads version x of the output of the berlind group finder and saves a list of galaxy 
    and group IDs. This output file is used to create useable group catalogues.

process_mockruns_4.0.py
    obsolete.  This was used to create a value added group catalogue for the berlind group 
    finder run over the mock.  This has been replaced by a script that does this for all 
    group finders.

process_all_mockruns.sh
    obsolete.  This ran the above script for all the mocks we have

bootstrap_mock_groups_2.0.py
    obsolete.  creates bootstrap group catalogues for the mock berlind group catalogues.

make_berlind_mock_bootstraps_4.0.sh
    obsolete.  runs the above script for all mocks.
    
all the code that makes useful group catalogues and bootstrap catalogues has been moved
to the galactic conformity project folder.