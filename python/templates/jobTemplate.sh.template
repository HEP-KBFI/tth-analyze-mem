#!/bin/bash -l

set -v

echo -n "Current time: "
date
echo -n "Hostname: "
hostname

echo -e "\n\nInitializing CMSSW run-time environment"

shopt -s expand_aliases
source /cvmfs/cms.cern.ch/cmsset_default.sh

echo -e "\n\nGoing to directory {{ cmsswSrcDir }} and calling cmsenv there"
cd {{ cmsswSrcDir }}
cmsenv

echo -e "\n\nCreating folder $(dirname {{ inFileNameScratch }})"
mkdir -p $(dirname {{ inFileNameScratch }})

echo -e "\n\nCreating folder $(dirname {{ outFileNameScratch }})"
mkdir -p $(dirname {{ outFileNameScratch }})

echo -e "\n\nCopying {{ inFileNameLocal }} to {{ inFileNameScratch }}"
cp {{ inFileNameLocal }} {{ inFileNameScratch }}

echo -e "\n\nRunning {{ execName }}"
echo "{{ execName }} {{ pythonCfg }}"
LINE=$(seq -s= 30|tr -d '[:digit:]')
echo "${LINE} BEGIN ${LINE}"
{{ execName }} {{ pythonCfg }}
echo "${LINE}= END =${LINE}"

echo -n -e "\n\nFinished at "
date

echo -e "\n\nCreating directory $(dirname {{ outFileNameLocal }})"
mkdir -p $(dirname {{ outFileNameLocal }})

echo -e "\n\nMoving {{ outFileNameScratch }} to {{ outFileNameLocal }}"
mv {{ outFileNameScratch }} {{ outFileNameLocal }}

echo -e "\n\nCleaning up"
rm -rf $(dirname {{ inFileNameScratch }})

echo -e "\n\nDone"

