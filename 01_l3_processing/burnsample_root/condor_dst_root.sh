#!/bin/bash

#eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`
#eval '/cvmfs/icecube.opensciencegrid.org/standard/icetray-start'
#eval cvmfs
#eval combo

# Define the output directory
OUTPUT_DIR="/data/user/gagrawal/unblind_root/"

TMPDIR=$(mktemp -d)
cp l3_i3_to_root_day.py $TMPDIR

cd $TMPDIR
python l3_i3_to_root_day.py $1 $2 $3 $4

if [ ! -d "$OUTPUT_DIR" ]; then
  mkdir -p "$OUTPUT_DIR"
  echo "Directory $OUTPUT_DIR created."
else
  echo "Directory $OUTPUT_DIR already exists."
fi

mv l3_data* $OUTPUT_DIR



