#!/bin/bash
#
#   Returns info about the current version of ARRM_V2
#
echo "ARRM_V2 version (from git)"        > ARRM_V2_version.txt      
git status | head -n 1                  >> ARRM_V2_version.txt      # current branch
git log | head -n 1                     >> ARRM_V2_version.txt      # commit hash
git log | head -n 3 | tail -n 1         >> ARRM_V2_version.txt      # timestamp of latest pull
git remote -v | head -n 1               >> ARRM_V2_version.txt      # URL of pull

mods=$(git status --porcelain)
if [[ $mods ]]
then
    echo  -e  "\nwith the following changes:\n" >> ARRM_V2_version.txt
    git status --porcelain >> ARRM_V2_version.txt
    exit 1
fi

exit 0

