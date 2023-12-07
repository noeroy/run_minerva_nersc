
#Configuration file to use
config_file="campaign_param.txt"

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup dunesw v09_45_00_00 -q e20:prof
setup edepsim v3_2_0 -q e20:prof

export WORKING_DIRECTORY=$(pwd)

export runnum=$1

#SETING UP THE CAMPAIGN PARAMETERS
CAMPAIGN=""
POT=""
MODE=""
FILELOCATION=""
OFFSETX=0
OFFSETY=0
OFFSETZ=0

if [[ -f "$config_file" ]]; then
  source "$config_file"
else
  echo "Warning: Configuration file not found."
fi


run=$( printf '%05d\n' $runnum)

export EDEPSIM_NAME="${CAMPAIGN}_${POT}_${MODE}.spill.${run}.EDEPSIM_SPILLS.root"
echo $runnum $run
export OUTPUT_DIRECTORY=$WORKING_DIRECTORY'/output_dir/run_'$run
export LOG_DIRECTORY=$OUTPUT_DIRECTORY'/logs'
export OPTION_DIRECTORY=$OUTPUT_DIRECTORY'/options'
export DST_DIRECTORY=$OUTPUT_DIRECTORY'/dst'
export FLATROOT_DIRECTORY=$OUTPUT_DIRECTORY'/flat_root'
export TMP_MINERVA=$OUTPUT_DIRECTORY'/minerva_root'
export FLATROOT_NAME="${CAMPAIGN}_${POT}_${MODE}.spill.${run}.EDEPSIM_SPILLS.flat.root"

which root
which python3

mkdir -p $OUTPUT_DIRECTORY
mkdir $LOG_DIRECTORY
mkdir $OPTION_DIRECTORY
mkdir $DST_DIRECTORY
mkdir $FLATROOT_DIRECTORY
mkdir $TMP_MINERVA


CONF_FILE="${OPTION_DIRECTORY}/conf_${run}.ini"
rm -f $CONF_FILE

echo '[data] '>> $CONF_FILE
echo "offset_x = ${OFFSETX} ">> $CONF_FILE
echo "offset_y = ${OFFSETY} ">> $CONF_FILE
echo "offset_z = ${OFFSETZ} ">> $CONF_FILE

echo "file_path = ${FILELOCATION} ">> $CONF_FILE
echo "file_name = ${EDEPSIM_NAME} ">> $CONF_FILE
echo "run_number = ${runnum} ">> $CONF_FILE
echo "output_dir = ${FLATROOT_DIRECTORY} ">> $CONF_FILE
echo "output_name = ${FLATROOT_NAME} ">> $CONF_FILE


export NSPILL=$(python3 convert_flat_edepsim.py --config $CONF_FILE 2>&1  | tail -n 1)

# Create option file 
DSTNAME=${CAMPAIGN}_${POT}_${MODE}.minerva.${run}.dst.root
NAME_ROOT_OUTPUT=${CAMPAIGN}_${POT}_${MODE}.minerva.${run}.IDODDigits.root
NAME_ROOT_HISTO=${CAMPAIGN}_${POT}_${MODE}.minerva.${run}.Histogam.root

FULLFILE=$OPTION_DIRECTORY'/full_'$run'.opts'

cp MC_full_minerva.opts $FULLFILE

sed -i "s/MAXEVT/${NSPILL}/g" $FULLFILE
sed -i "s#TMP_MINERVA#${TMP_MINERVA}#g" $FULLFILE
sed -i "s#FLATROOT_DIRECTORY#${FLATROOT_DIRECTORY}#g" $FULLFILE
sed -i "s#FLATROOT_NAME#${FLATROOT_NAME}#g" $FULLFILE
sed -i "s#DST_DIRECTORY#${DST_DIRECTORY}#g" $FULLFILE

sed -i "s/DSTNAME/${DSTNAME}/g" $FULLFILE
sed -i "s#NAME_ROOT_OUTPUT#${NAME_ROOT_OUTPUT}#g" $FULLFILE
sed -i "s#NAME_ROOT_HISTO#${NAME_ROOT_HISTO}#g" $FULLFILE

