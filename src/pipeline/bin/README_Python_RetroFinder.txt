Running the Python RetroFinder pipeline:

NOTE: Provide full path to config file as some scripts rely on that as can 
not assume it is in the same directory from which the script is being run.
Make sure the config file has the correct database, data version and date.
Requires Mark's pycbio directory of Python libraries to be symlinked to the
retroFinder/trunk/src/pipeline/lib directory
cd <path>/retroFinder/trunk/src/pipeline/lib
ln -s /hive/groups/gencode/local/pycbio .

Scripts to run:
getSequenceData <path>/configFile.cfg
prepareSeqsForAlignment <path>/configFile.cfg
alignMrnaSeqs <path>/configFile.cfg
filterAndChainMrnaAlignments <path>/configFile.cfg
getAnnotForRetroPred <path>/configFile.cfg
doRetrogenePrediction <path>/configFile.cfg
