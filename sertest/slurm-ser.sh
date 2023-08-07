#!/bin/bash

#SBATCH --partition=priority
#SBATCH --constraint='skylake'
#SBATCH --qos=pol02003sky
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=snakeser
#SBATCH -o ser-%j.out
#SBATCH -e ser-%j.err

LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/lib"
export TIMEFORMAT="user-seconds %3U"
cd /home/pol02003/dub/sertest
time ./dubser
