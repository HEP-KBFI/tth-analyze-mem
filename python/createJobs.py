import os, codecs, sys, logging, getpass, stat
from tthAnalysis.tthMEM.jobTemplates import getNofEntries, \
  createPythonCfg, createBashCfg, createSbatch, createMakefile, createPythonROCcfg

def createJobs(samples, channel, version, lepton_selections, central_or_shifts, execName,
               treeName, integrationMode, maxObjFunctionCalls, nofIntegrationsPerJob,
               lhRatioBranchName, rocLegendPosition):
  '''
  TODO: - remove unnecessary complexity in the paths (currently both the file name and dirname
          contain the same information)
        - make the file names of roc curve plots channel and version specific
  '''

  if os.environ.get('CMSSW_BASE') is None:
    logging.error("Variable CMSSW_BASE unset! " + \
                  "Execute the script only if you have a working CMSSW space set up")
    sys.exit(1)

  baseDirPattern = os.path.join("/%s", getpass.getuser(), "ttHAnalysis", version)
  baseDir = baseDirPattern % "home"
  scratchDir = baseDirPattern % "scratch"
  cmsswSrcDir = os.path.join(os.environ.get('CMSSW_BASE'), "src")
  rocDir = os.path.join(baseDir, "mem", "roc")
  rocPlotDir = os.path.join(rocDir, "plots")
  rocCmd = "memROC"

  sbatchBashFiles, sbatchLogFiles = [], []
  outFileNameLocalArray = {}
  rocLabels, rocOutFileNames = [], []
  inputSignalFile = ""
  inputBkgFiles = []

  for sampleKey, sampleValue in samples.iteritems():
    if not sampleValue["use_it"]: continue

    if sampleValue["sample_category"] == "signal":
      rocLabels.insert(0, sampleValue["process_name_specific"])
    else:
      rocLabels.append(sampleValue["process_name_specific"])

    for lepton_selection in lepton_selections:
      for central_or_shift in central_or_shifts:

        filePattern = os.path.join(
          "%s", "%s", channel, lepton_selection, sampleValue["process_name_specific"],
          "_".join(["%s", channel, sampleValue["process_name_specific"],
                    lepton_selection, central_or_shift])
        )

        fileNameLocal = filePattern % (baseDir, "output_root", "out") + ".root"
        if not os.path.exists(fileNameLocal) or not os.path.isfile(fileNameLocal):
          logging.warning("File %s does not exists! Skipping" % fileNameLocal)
          continue
        else:
          logging.info("Create configuration files for %s" % fileNameLocal)

        fileNameScratch = filePattern % (scratchDir, os.path.join("output_root", "%d"), "out") + ".root"
        outFileNameLocal = filePattern % (baseDir, os.path.join("mem", "output"), "out") + "_%d.root"
        outFileNameScratch = filePattern % (scratchDir, os.path.join("mem", "output"), "mem") + "_%d.root"
        jobCfgFile = filePattern % (baseDir, os.path.join("mem", "cfg"), "cfg") + "_%d.py"
        jobBashFile = filePattern % (baseDir, os.path.join("mem", "cfg"), "cfg") + "_%d.sh"
        logFile = filePattern % (baseDir, os.path.join("mem", "log"), "log") + "_%d.txt"

        outFileNameLocalResult = filePattern % (baseDir, os.path.join("mem", "output"), "mem") + ".root"
        outFileNameLocalArray[outFileNameLocalResult] = []
        
        if inputSignalFile == "" and sampleValue["sample_category"] == "signal":
          inputSignalFile = outFileNameLocalResult
        elif outFileNameLocalResult not in inputBkgFiles:
          inputBkgFiles.append(outFileNameLocalResult)

        jobDir = os.path.dirname(jobCfgFile)
        logDir = os.path.dirname(logFile)
        for d in [jobDir, logDir]:
          if not os.path.exists(d) or not os.path.isdir(d):
            os.makedirs(d)

        nofEntries = getNofEntries(fileNameLocal, treeName)
        if nofEntries == 0:
          logging.warning("Tree %s in file %s contains no events. Skipping" % (treeName, fileNameLocal))
          sys.exit(1)

        entryStartingPoints = range(0, nofEntries, nofIntegrationsPerJob)
        nofJobs = len(entryStartingPoints)
        for i in range(nofJobs):
          fileNameScratch_i = fileNameScratch % i
          outFileNameLocal_i = outFileNameLocal % i
          outFileNameScratch_i = outFileNameScratch % i
          jobCfgFile_i = jobCfgFile % i
          jobBashFile_i = jobBashFile % i
          logFile_i = logFile % i

          startingPoint = entryStartingPoints[i]
          nofEventsToProcess = nofIntegrationsPerJob if i != nofJobs - 1 else -1

          pythonCfg = createPythonCfg(
            fileNameScratch_i, nofEventsToProcess, outFileNameScratch_i, treeName,
            integrationMode, maxObjFunctionCalls, startingPoint
          )
          bashCfg = createBashCfg(
            fileNameLocal, outFileNameLocal_i, fileNameScratch_i, outFileNameScratch_i,
            execName, jobCfgFile_i, cmsswSrcDir
          )
          with codecs.open(jobCfgFile_i, 'w', 'utf-8') as f: f.write(pythonCfg)
          with codecs.open(jobBashFile_i, 'w', 'utf-8') as f: f.write(bashCfg)
          st = os.stat(jobBashFile_i)
          os.chmod(jobBashFile_i, st.st_mode | stat.S_IEXEC)

          sbatchBashFiles.append(jobBashFile_i)
          sbatchLogFiles.append(logFile_i)
          outFileNameLocalArray[outFileNameLocalResult].append(outFileNameLocal_i)

  for i in range(1, len(rocLabels)):
    rocOutFileNames.append(os.path.join(rocPlotDir, "roc_%s_%s.pdf" % (rocLabels[0], rocLabels[i])))
  if len(rocLabels) > 2:
    rocOutFileNames.append(os.path.join(rocPlotDir, "roc.pdf"))

  sbatchFile = os.path.join(baseDir, "mem", "sbatchMEM_%s.py" % channel)
  sbatchContents = createSbatch(sbatchBashFiles, sbatchLogFiles)
  with codecs.open(sbatchFile, 'w', 'utf-8') as f: f.write(sbatchContents)
  st = os.stat(sbatchFile)
  os.chmod(sbatchFile, st.st_mode | stat.S_IEXEC)

  rocCfgFile = os.path.join(rocDir, "cfg.py")
  if not os.path.exists(rocDir):
    os.makedirs(rocDir)
  rocCfgContents = createPythonROCcfg(inputSignalFile, inputBkgFiles, rocPlotDir, treeName,
                                      lhRatioBranchName, rocLabels, rocLegendPosition)
  with codecs.open(rocCfgFile, 'w', 'utf-8') as f: f.write(rocCfgContents)

  makeFile = os.path.join(baseDir, "_".join(["Makefile", channel, "mem"]))
  makeFileContents = createMakefile(sbatchFile, outFileNameLocalArray, scratchDir,
                                    rocOutFileNames, rocCmd, rocCfgFile)
  with codecs.open(makeFile, 'w', 'utf-8') as f: f.write(str(makeFileContents))

  logging.info("Run:\tmake -f %s -j 4" % makeFile)
  logging.info("Done")
