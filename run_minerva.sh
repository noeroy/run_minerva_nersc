export WORKING_DIRECTORY=$(pwd)

export runnum=$1

run=$( printf '%05d\n' $runnum)

echo $runnum $run
export OUTPUT_DIRECTORY=$WORKING_DIRECTORY'/output_dir/run_'$run
export LOG_DIRECTORY=$OUTPUT_DIRECTORY'/logs'
export OPTION_DIRECTORY=$OUTPUT_DIRECTORY'/options'
export DST_DIRECTORY=$OUTPUT_DIRECTORY'/dst'
export FLATROOT_DIRECTORY=$OUTPUT_DIRECTORY'/flat_root'
export TMP_MINERVA=$OUTPUT_DIRECTORY'/minerva_root'

echo "SETUP MINERVA"

#for runs before pico 4.1
#source setup_minerva_pico4.1.sh

#Updated geometry
source setup_minerva.sh

cd $WORKING_DIRECTORY

echo $run
echo "FULL JOB"
echo $(date)
SystemTestsApp.exe $OPTION_DIRECTORY'/full_'$run'.opts' > $LOG_DIRECTORY'/log_full_'$run'.logs'

echo " "
echo "done."
echo $(date)

