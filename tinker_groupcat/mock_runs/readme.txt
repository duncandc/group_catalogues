what do all these scripts do?

read_tinker_mockruns.py
    reads the output of the yang group finder and saves a list of galaxy and group IDs.
    This output file is used to create useable group catalogues.
    
read_all_tinier_mockruns.sh
    runs the script above for all the mocks we have.
    
process_tinier_mockruns.py
    obsolete.  This was used to create a value added group catalogue for the yang group 
    finder run over the mock.  This has been replaced by a script that does this for all 
    group finders.

process_all_tinker_mockruns.sh
    obsolete.  This ran the above script for all the mocks we have
    
bootstrap_mock_groups.py
    creates bootstrap group catalogues for the mock yang group catalogues.

make_tinker_mock_bootstraps.sh
    runs the above script for all mocks.