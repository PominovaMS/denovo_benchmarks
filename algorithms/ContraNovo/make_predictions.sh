#!/bin/bash

# Activate conda environment
conda init
source /opt/conda/etc/profile.d/conda.sh
conda activate ContraNovo

# Get dataset property tags
DSET_TAGS=$(python /algo/base/dataset_tags_parser.py --dataset "$@")
while IFS='=' read -r key value; do
    export "$key"="$value"
done <<< "$DSET_TAGS"

# Choose the device  
if command -v nvidia-smi &> /dev/null
then
    nvidia-smi > /dev/null 2>&1
    if [ $? -eq 0 ]; then
        device=0 # GPU:0
    else
        device=-1 # CPU
    fi
else
    device=-1 # CPU
fi
echo "Using device $device."

# Use tag variables to specify de novo algorithm
# for the particular dataset properties
cd /algo/ContraNovo
echo "Processing mgf files:"
python -m ContraNovo.ContraNovo --mode=denovo --peak_path="/algo/$@"/*.mgf --model=./ContraNovo/ContraNovo.ckpt

# Convert predictions to the general output format
cd /algo
mv ./ContraNovo/outputs.mztab ./
python ./output_mapper.py --output_path=outputs.mztab --input_dir="$@"