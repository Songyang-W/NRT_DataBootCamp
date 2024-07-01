#!/bin/bash -l

#$ -P jchenlab
#$ -j y               # Merge the error and output streams into a single file
#$ -pe omp 16
#$ -l num_proc=16
#$ -m bea

#$ -t 1
 
ZDRIVE=/net/claustrum/mnt/data
YDRIVE=/net/claustrum/mnt/data1
XDRIVE=/net/claustrum2/mnt/data
# program name or command and its options and arguments
module load matlab
cd /usr3/graduate/dglee3/2P/
directory='/net/claustrum/mnt/data/Projects/Perirhinal/Animals/'
anm=$1
sess=$2

mkdir -p /scratch/chen_step2_prh/"$anm"-"$sess"/PreProcess/A0_Ch0
mkdir -p /scratch/chen_step2_prh/"$anm"-"$sess"/PreProcess/A1_Ch0

matlab -nodisplay -r "poolobj=gcp('nocreate'); delete(poolobj); pr_pipeline_Step2_SCC('$anm','$sess'); poolobj=gcp('nocreate'); delete(poolobj); exit"

chmod -R g+w "/net/claustrum2/mnt/data/Projects/Perirhinal/Animals/$anm/2P/$anm-$sess/PreProcess/A0_Ch0"
chmod -R g+w "/net/claustrum2/mnt/data/Projects/Perirhinal/Animals/$anm/2P/$anm-$sess/PreProcess/A1_Ch0"
chmod g+w "/net/claustrum2/mnt/data/Projects/Perirhinal/Animals/$anm/$anm-$sess".mat
mkdir -p /scratch/chen_step2_prh/"$anm"-"$sess"/PreProcess/"$area" && chmod g+w -R /scratch/chen_step2_prh