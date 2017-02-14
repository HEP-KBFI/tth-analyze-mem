import os, codecs, sys, logging, getpass, stat
from tthAnalysis.tthMEM.jobTemplates import getNofEntries, \
  createPythonCfg, createBashCfg, createSbatch, createMakefile, createPythonROCcfg

class JobCreator:
  def __init__(self, samples, channel, year, version, memBaseDir, central_or_shifts, charge_selections,
               lepton_selections, execName, treeName, rleSelectionFile, integrationMode,
               maxObjFunctionCalls, nofIntegrationsPerJob, lhRatioBranchName, debugPlots, forceGenLevel,
               higgsWidth, clampVariables, markovChainParams, comment, priority = 'main'):
    self.samples               = samples
    self.channel               = channel
    self.year                  = year
    self.version               = version
    self.memBaseDir            = memBaseDir
    self.central_or_shifts     = central_or_shifts
    self.charge_selections     = charge_selections
    self.lepton_selections     = lepton_selections
    self.execName              = execName
    self.treeName              = treeName
    self.rleSelectionFile      = rleSelectionFile
    self.integrationMode       = integrationMode
    self.maxObjFunctionCalls   = maxObjFunctionCalls
    self.nofIntegrationsPerJob = nofIntegrationsPerJob
    self.lhRatioBranchName     = lhRatioBranchName
    self.debugPlots            = debugPlots
    self.forceGenLevel         = forceGenLevel
    self.higgsWidth            = higgsWidth
    self.clampVariables        = clampVariables
    self.markovChainParams     = markovChainParams
    self.comment               = comment
    self.priority              = priority

    self.baseDirPattern = os.path.join(getpass.getuser(), "ttHAnalysis", self.year, self.version)
    self.baseDir        = os.path.join("/home", self.baseDirPattern)
    self.cmsswSrcDir    = os.path.join(os.environ.get('CMSSW_BASE'), "src")

    self.scratchDir           = os.path.join("/scratch", self.baseDirPattern)
    self.scratchTempOutputDir = os.path.join(self.scratchDir, "temp_output", "%d")
    self.scratchOutputDir     = os.path.join(self.scratchDir, "mem_output")

    self.memDir         = os.path.join(self.baseDir, self.memBaseDir)
    self.memCfgDir      = os.path.join(self.memDir, "cfg")
    self.memOutputDir   = os.path.join(self.memDir, "output")
    self.memLogDir      = os.path.join(self.memDir, "log")
    self.memCommentFile = os.path.join(self.memDir, "comment.txt")

    self.rocDir     = os.path.join(self.memDir, "roc")
    self.rocPlotDir = os.path.join(self.rocDir, "plots")
    self.rocCSVDir  = os.path.join(self.rocDir, "csvs")
    self.rocCmd     = "memROC"
    self.rocCfgFile = os.path.join(self.rocDir, "cfg.sh")

    self.sbatchFile = os.path.join(self.memDir, "sbatchMEM_%s.py" % self.channel)
    self.makeFile   = os.path.join(self.memDir, "_".join(["Makefile", self.channel, "mem"]))

  @staticmethod
  def isMEMPathAvailable(year, version, memBaseDir):
    return not os.path.isdir(os.path.join("/home", getpass.getuser(), "ttHAnalysis", year, version, memBaseDir))

  def createJobs(self):

    if os.environ.get('CMSSW_BASE') is None:
      logging.error("Variable CMSSW_BASE unset! " + \
                    "Execute the script only if you have a working CMSSW space set up")
      sys.exit(1)

    sbatchBashFiles, sbatchLogFiles, sbatchOutFileNameLocalFiles = [], [], []
    outFileNameLocalArray = {}
    rocLabels, rocOutFileNames, rocCSVOutFileNames = [], [], []
    inputSignalFile = ""
    inputBkgFiles = {}

    for sampleValue in self.samples.values():
      if not sampleValue["use_it"]: continue

      process_name = sampleValue["process_name"]
      if sampleValue["sample_category"] == "signal":
        rocLabels.insert(0, process_name)
      elif process_name not in rocLabels:
        rocLabels.append(process_name)

      isMC = (sampleValue["type"] == "mc")
      is2016 = (self.year == "2016")

      for lepton_selection in self.lepton_selections:
        for charge_selection in self.charge_selections:
          for central_or_shift in self.central_or_shifts:
            for sample_name in sampleValue["process_name_specific"][self.year]:

              filePatternGen = "_".join(["%s", self.channel, process_name, lepton_selection, charge_selection, central_or_shift])
              filePattern    = "_".join(["%s", self.channel, sample_name,  lepton_selection, charge_selection, central_or_shift])
              fileNameLocal  = os.path.join(
                self.baseDir, "output_root", self.channel, "_".join([lepton_selection, charge_selection]), sample_name, filePattern % "out"
              ) + ".root"

              if not os.path.exists(fileNameLocal) or not os.path.isfile(fileNameLocal):
                logging.warning("File %s does not exists! Skipping" % fileNameLocal)
                continue
              else:
                logging.info("Create configuration files for %s" % fileNameLocal)

              fileNameScratch           = os.path.join(self.scratchTempOutputDir, filePattern    % "out") + ".root"
              outFileNameLocal          = os.path.join(self.memOutputDir,         filePattern    % "out") + "_%d.root"
              outFileNameScratch        = os.path.join(self.scratchOutputDir,     filePattern    % "mem") + "_%d.root"
              jobCfgFile                = os.path.join(self.memCfgDir,            filePattern    % "cfg") + "_%d.py"
              jobBashFile               = os.path.join(self.memCfgDir,            filePattern    % "cfg") + "_%d.sh"
              logFile                   = os.path.join(self.memLogDir,            filePattern    % "log") + "_%d.txt"
              outFileNameLocalResult    = os.path.join(self.memOutputDir,         filePattern    % "out") + ".root"
              outFileNameLocalResultGen = os.path.join(self.memOutputDir,         filePatternGen % "out") + ".root"
              outFileNameLocalArray[outFileNameLocalResult] = []

              if inputSignalFile == "" and sampleValue["sample_category"] == "signal":
                inputSignalFile = outFileNameLocalResult
              elif sampleValue["sample_category"] != "signal":
                if outFileNameLocalResultGen not in inputBkgFiles:
                  inputBkgFiles[outFileNameLocalResultGen] = []
                if outFileNameLocalResult not in inputBkgFiles[outFileNameLocalResultGen]:
                  inputBkgFiles[outFileNameLocalResultGen].append(outFileNameLocalResult)

              jobDir = os.path.dirname(jobCfgFile)
              logDir = os.path.dirname(logFile)
              for d in [jobDir, logDir]:
                if not os.path.exists(d) or not os.path.isdir(d):
                  os.makedirs(d)

              nofEntries = getNofEntries(fileNameLocal, self.treeName)
              if nofEntries == 0:
                logging.warning("Tree %s in file %s contains no events. Skipping" % (self.treeName, fileNameLocal))
                sys.exit(1)

              entryStartingPoints = range(0, nofEntries, self.nofIntegrationsPerJob)
              nofJobs = len(entryStartingPoints)
              for i in range(nofJobs):
                fileNameScratch_i    = fileNameScratch % i
                outFileNameLocal_i   = outFileNameLocal % i
                outFileNameScratch_i = outFileNameScratch % i
                jobCfgFile_i         = jobCfgFile % i
                jobBashFile_i        = jobBashFile % i
                logFile_i            = logFile % i

                startingPoint = entryStartingPoints[i]
                nofEventsToProcess = self.nofIntegrationsPerJob if i != nofJobs - 1 else -1

                pythonCfg = createPythonCfg(
                  isMC, is2016, fileNameScratch_i, nofEventsToProcess, outFileNameScratch_i, self.treeName,
                  self.integrationMode, self.maxObjFunctionCalls, startingPoint, self.debugPlots,
                  self.forceGenLevel, self.higgsWidth, self.clampVariables, self.markovChainParams, self.rleSelectionFile
                )
                bashCfg = createBashCfg(
                  fileNameLocal, outFileNameLocal_i, fileNameScratch_i, outFileNameScratch_i,
                  self.execName, jobCfgFile_i, self.cmsswSrcDir
                )
                with codecs.open(jobCfgFile_i, 'w', 'utf-8') as f: f.write(pythonCfg)
                with codecs.open(jobBashFile_i, 'w', 'utf-8') as f: f.write(bashCfg)
                st = os.stat(jobBashFile_i)
                os.chmod(jobBashFile_i, st.st_mode | stat.S_IEXEC)

                sbatchBashFiles.append(jobBashFile_i)
                sbatchLogFiles.append(logFile_i)
                sbatchOutFileNameLocalFiles.append(outFileNameLocal_i)
                outFileNameLocalArray[outFileNameLocalResult].append(outFileNameLocal_i)

    if not outFileNameLocalArray:
      sys.exit(0)

    for i in range(1, len(rocLabels)):
      rocOutFileNames.append(os.path.join(self.rocPlotDir, "roc_%s_%s.pdf" % (rocLabels[0], rocLabels[i])))
    if len(rocLabels) > 2:
      rocOutFileNames.append(os.path.join(self.rocPlotDir, "roc.pdf"))

    sbatchContents = createSbatch(sbatchBashFiles, sbatchLogFiles, sbatchOutFileNameLocalFiles, self.priority)
    with codecs.open(self.sbatchFile, 'w', 'utf-8') as f:
      f.write(sbatchContents)
    st = os.stat(self.sbatchFile)
    os.chmod(self.sbatchFile, st.st_mode | stat.S_IEXEC)


    if not os.path.exists(self.rocDir):
      os.makedirs(self.rocDir)
    allInputFiles = [inputSignalFile] + inputBkgFiles.keys()
    classLabels   = ['sig'] + (['bkg'] * len(inputBkgFiles))
    rocCfgContents = createPythonROCcfg(
      allInputFiles, classLabels, self.memDir, self.rocCSVDir, rocLabels[0], rocLabels[1:],
      rocOutFileNames, self.treeName, self.lhRatioBranchName
    )
    with codecs.open(self.rocCfgFile, 'w', 'utf-8') as f:
      f.write(rocCfgContents)
    st = os.stat(self.rocCfgFile)
    os.chmod(self.rocCfgFile, st.st_mode | stat.S_IEXEC)

    makeFileContents = createMakefile(self.sbatchFile, outFileNameLocalArray, self.scratchDir,
                                      rocOutFileNames, self.rocCfgFile, inputBkgFiles, inputSignalFile)
    with codecs.open(self.makeFile, 'w', 'utf-8') as f:
      f.write(str(makeFileContents))

    if self.comment:
      with codecs.open(self.memCommentFile, 'w', 'utf-8') as f:
        f.write(self.comment)
