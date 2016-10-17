import subprocess, jinja2

rootEventCounterTemplate = r"""TFile * f = new TFile(\"{{ rootFile }}\", \"read\"); \
TTree * t = static_cast<TTree *>(f -> Get(\"{{ treeName }}\")); \
std::cout << t -> GetEntries();"""

pythonCfgTemplate = r"""import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
  fileNames = cms.vstring('{{ inFileName }}'),
  maxEvents = cms.int32({{ maxEvents }}),
  outputEvery = cms.uint32(50)
)

process.fwliteOutput = cms.PSet(
  fileName = cms.string('{{ outFileName }}')
)

process.logging = cms.PSet(
  logLevel = cms.string("info"),
  enableTimeStamp = cms.bool(True),
  enableLogging = cms.bool(True)
)

process.tthMEM = cms.PSet(
  isMC                = cms.bool({{ isMC }}),
  is2016              = cms.bool({{ is2016 }}),
  treeName            = cms.string("{{ treeName }}"),
  rleSelectionFile    = cms.string("{{ rleSelectionFile  }}"),
  pdfName             = cms.string("MSTW2008lo68cl"),
  madgraphFileName    = cms.string("tthAnalysis/tthMEM/data/param_card.dat"),
  integrationMode     = cms.string("{{ integrationMode }}"),
  maxObjFunctionCalls = cms.uint32({{ maxObjFunctionCalls }}),
  startingFromEntry   = cms.int64({{startingFromEntry }}),
  debugPlots          = cms.uint32({{ debugPlots }}),
  forceGenLevel       = cms.bool({{ forceGenLevel }}),
  higgsWidth          = cms.double({{ higgsWidth }}),
  clampVariables      = cms.VPSet({% for it in clampVariables %}
    cms.PSet( var = cms.string("{{ it.0 }}"), useGen = cms.bool({{ it.1 }}), useCfg = cms.bool({{ it.2 }}), val = cms.double({{ it.3 }})),{% endfor %}
  ),
  markovChainParams   = cms.PSet(
    mode                = cms.string("{{ markovChainParams["mode"] }}"),
    nofBatches          = cms.uint32({{ markovChainParams["nofBatches"] }}),
    nofChains           = cms.uint32({{ markovChainParams["nofChains"] }}),
    maxCallsStartingPos = cms.uint32({{ markovChainParams["maxCallsStartingPos"] }}),
    epsilon0            = cms.double({{ markovChainParams["epsilon0"] }}),
    T0                  = cms.double({{ markovChainParams["T0"] }}),
    nu                  = cms.double({{ markovChainParams["nu"] }})
  )
)

"""

pythonROCcfgTemplate = r"""import FWCore.ParameterSet.Config as cms

roc = cms.PSet()

roc.tthMEM = cms.PSet(
  signalFile = cms.string('{{ signalFile }}'),
  bkgFiles = cms.vstring('{{ bkgFiles|join('\', \'') }}'),
  outFolder = cms.string('{{ outFolder }}'),
  treeName = cms.string('{{ treeName }}'),
  branchName = cms.string('{{ branchName }}'),
  labels = cms.vstring('{{ labels|join('\', \'') }}'),
  legendPosition = cms.vdouble({{ legendPosition|join(', ') }})
)

"""

jobTemplate = r"""#!/bin/bash

echo -n "Current time: "
date
echo -n "Hostname: "
hostname

echo -e "\n\nInitializing CMSSW run-time environment"

echo "shopt -s expand_aliases"
shopt -s expand_aliases

echo "source /cvmfs/cms.cern.ch/cmsset_default.sh"
source /cvmfs/cms.cern.ch/cmsset_default.sh

echo -e "\n\nGoing to directory {{ cmsswSrcDir }} and calling cmsenv there"
echo "cd {{ cmsswSrcDir }}"
cd {{ cmsswSrcDir }}
cmsenv

echo -e "\n\nCreating folder $(dirname {{ inFileNameScratch }})"
echo "mkdir -p $(dirname {{ inFileNameScratch }})"
mkdir -p $(dirname {{ inFileNameScratch }})

echo -e "\n\nCreating folder $(dirname {{ outFileNameScratch }})"
echo "mkdir -p $(dirname {{ outFileNameScratch }})"
mkdir -p $(dirname {{ outFileNameScratch }})

echo -e "\n\nCopying {{ inFileNameLocal }} to {{ inFileNameScratch }}"
echo "cp {{ inFileNameLocal }} {{ inFileNameScratch }}"
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
echo "mkdir -p $(dirname {{ outFileNameLocal }})"
mkdir -p $(dirname {{ outFileNameLocal }})

echo -e "\n\nMoving {{ outFileNameScratch }} to {{ outFileNameLocal }}"

echo "mkdir -p $(dirname {{ outFileNameLocal }})"
mkdir -p $(dirname {{ outFileNameLocal }})

echo "mv {{ outFileNameScratch }} {{ outFileNameLocal }}"
mv {{ outFileNameScratch }} {{ outFileNameLocal }}

echo -e "\n\nCleaning up"
echo "rm -rf $(dirname {{ inFileNameScratch }})"
rm -rf $(dirname {{ inFileNameScratch }})

echo -e "\n\nDone"

"""

sbatchTemplate = r"""#!/usr/bin/env python
import subprocess, sys, getpass, time

if __name__ == '__main__':
  commands, sbatchIDs = [], []
  {% for bashScript, logFile in zippedScriptLog %}
  commands.append('sbatch --output={{ logFile }} {{ bashScript }}') {% endfor %}
  for command in commands:
    submitJobProcess = subprocess.Popen(command,
                                        stdout = subprocess.PIPE,
                                        stderr = subprocess.PIPE,
                                        shell = True)
    submitStdout, submitStderr = submitJobProcess.communicate()
    sbatchIDs.append(submitStdout.rstrip('\n').split()[-1])

  squeueCommand = "squeue -u %s | tail -n+2 | awk '{print $1}'" % getpass.getuser()
  while True:
    squeueProcess = subprocess.Popen(squeueCommand,
                                     stdout = subprocess.PIPE,
                                     stderr = subprocess.PIPE,
                                     shell = True)
    squeueStdout, squeueStderr = squeueProcess.communicate()
    nofJobsStillRunning = 0
    for line in squeueStdout.rstrip('\n').split('\n'):
      if line in sbatchIDs: nofJobsStillRunning += 1
    if nofJobsStillRunning == 0: break
    else:                        time.sleep(30)

"""

makefileTemplate = """
.DEFAULT_GOAL := all
SHELL := /bin/bash

run:
\t{{ waitingScript }}
{% for outFileNameLocalResult, outFileNameLocals in outFileNameLocalArray.iteritems() %}
{% for outFileNameLocal in outFileNameLocals %}
{{ outFileNameLocal }}: run
\t:
{% endfor %}
{{ outFileNameLocalResult }}: {{ outFileNameLocals|join(' ') }}
\thadd -f {{ outFileNameLocalResult }} {{ outFileNameLocals|join(' ') }} {% endfor %}

{% for bkgParent, bkgArray in inputBkgFiles.iteritems() %}{% if bkgArray|length > 1 %}
{{ bkgParent }}: {{ bkgArray|join(' ') }}
\thadd -f {{ bkgParent }} {{ bkgArray|join(' ') }}{% endif %}{% endfor %}

{% for outFileName in rocOutFileNames %}
{{ outFileName }}: {{ inputBkgFiles.keys()|join(' ') }} {{ inputSignalFile }}
\t{{ rocCmd }} {{ rocCfg }}
{% endfor %}

all: {{ rocOutFileNames|join(' ') }}

.PHONY: clean

clean:{% for outFileNameLocalResult, outFileNameLocals in outFileNameLocalArray.iteritems() %}
\trm -f {{ outFileNameLocalResult }} {% for outFileNameLocal in outFileNameLocals %}
\trm -f {{ outFileNameLocal }} {% endfor %} {% endfor %} {% for outFileName in rocOutFileNames %}
\trm -f {{ outFileName }} {% endfor %}
\trm -rf {{ scratchDir }}

"""

def getNofEntries(fileName, treeName):
  rootEventCounterCode = jinja2.Template(rootEventCounterTemplate).render(
                            rootFile = fileName,
                            treeName = treeName)
  rootEventCounterCommand = "echo \"%s\" | root -b -l" % rootEventCounterCode
  rootEventCounterProcess = subprocess.Popen(rootEventCounterCommand,
                                             stdout = subprocess.PIPE,
                                             stderr = subprocess.PIPE,
                                             shell = True)
  counterStdout, counterStderr = rootEventCounterProcess.communicate()
  if counterStderr != "": raise ValueError(counterStderr)
  return int(counterStdout)

def createPythonCfg(isMC, is2016, inFileName, maxEvents, outFileName, treeName,
                    integrationMode, maxObjFunctionCalls, startingFromEntry,
                    debugPlots, higgsWidth, clampVariables, markovChainParams,
                    rleSelectionFile, forceGenLevel):
  return jinja2.Template(pythonCfgTemplate).render(
    isMC = isMC,
    is2016 = is2016,
    inFileName = inFileName,
    maxEvents = maxEvents,
    outFileName = outFileName,
    treeName = treeName,
    rleSelectionFile = rleSelectionFile,
    integrationMode = integrationMode,
    maxObjFunctionCalls = maxObjFunctionCalls,
    startingFromEntry = startingFromEntry,
    debugPlots = debugPlots,
    forceGenLevel = forceGenLevel,
    higgsWidth = higgsWidth,
    clampVariables = clampVariables,
    markovChainParams = markovChainParams)

def createPythonROCcfg(signalFile, bkgFiles, outFolder, treeName,
                       branchName, labels, legendPosition):
  return jinja2.Template(pythonROCcfgTemplate).render(
    signalFile = signalFile,
    bkgFiles = bkgFiles,
    outFolder = outFolder,
    treeName = treeName,
    branchName = branchName,
    labels = labels,
    legendPosition = legendPosition)

def createBashCfg(inFileNameLocal, outFileNameLocal, inFileNameScratch,
                  outFileNameScratch, execName, pythonCfg, cmsswSrcDir):
  return jinja2.Template(jobTemplate).render(
    inFileNameLocal = inFileNameLocal,
    outFileNameLocal = outFileNameLocal,
    inFileNameScratch = inFileNameScratch,
    outFileNameScratch = outFileNameScratch,
    execName = execName,
    pythonCfg = pythonCfg,
    cmsswSrcDir = cmsswSrcDir)

def createSbatch(bashScript, logFile):
  return jinja2.Template(sbatchTemplate).render(
    zippedScriptLog = zip(bashScript, logFile))

def createMakefile(waitingScript, outFileNameLocalArray, scratchDir,
                   rocOutFileNames, rocCmd, rocCfg, inputBkgFiles, inputSignalFile):
  return jinja2.Template(makefileTemplate).render(
    waitingScript = waitingScript,
    outFileNameLocalArray = outFileNameLocalArray,
    scratchDir = scratchDir,
    rocOutFileNames = rocOutFileNames,
    rocCmd = rocCmd,
    rocCfg = rocCfg,
    inputBkgFiles = inputBkgFiles,
    inputSignalFile = inputSignalFile)
