import subprocess, jinja2, os, ROOT

templatesDir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'templates')

rootEventCounterTemplate = open(os.path.join(templatesDir, 'rootEventCounterTemplate.cc.template')).read()
pythonCfgTemplate        = open(os.path.join(templatesDir, 'pythonCfgTemplate.py.template')).read()
pythonROCcfgTemplate     = open(os.path.join(templatesDir, 'pythonROCcfgTemplate.py.template')).read()
jobTemplate              = open(os.path.join(templatesDir, 'jobTemplate.sh.template')).read()
sbatchTemplate           = open(os.path.join(templatesDir, 'sbatchTemplate.py.template')).read()
makefileTemplate         = open(os.path.join(templatesDir, 'Makefile.template')).read()

def getNofEntries(fileName, treeName, usePyROOT = True):
  if usePyROOT:
    rootFile = ROOT.TFile(fileName, "read")
    tree = rootFile.Get(treeName)
    return int(tree.GetEntries())
  else:
    rootEventCounterCode = jinja2.Template(rootEventCounterTemplate).render(
                              rootFile = fileName,
                              treeName = treeName)
    rootEventCounterCommand = "echo \"%s\" | root -b -l" % rootEventCounterCode
    rootEventCounterProcess = subprocess.Popen(rootEventCounterCommand,
                                               stdout = subprocess.PIPE,
                                               stderr = subprocess.PIPE,
                                               shell  = True)
    counterStdout, counterStderr = rootEventCounterProcess.communicate()
    if counterStderr != "": raise ValueError(counterStderr)
    return int(counterStdout)

def createPythonCfg(isMC, is2016, inFileName, maxEvents, outFileName, treeName,
                    integrationMode, maxObjFunctionCalls, startingFromEntry,
                    debugPlots, forceGenLevel, higgsWidth, clampVariables,
                    markovChainParams, rleSelectionFile):
  return jinja2.Template(pythonCfgTemplate).render(
    isMC                = isMC,
    is2016              = is2016,
    inFileName          = inFileName,
    maxEvents           = maxEvents,
    outFileName         = outFileName,
    treeName            = treeName,
    rleSelectionFile    = rleSelectionFile,
    integrationMode     = integrationMode,
    maxObjFunctionCalls = maxObjFunctionCalls,
    startingFromEntry   = startingFromEntry,
    debugPlots          = debugPlots,
    forceGenLevel       = forceGenLevel,
    higgsWidth          = higgsWidth,
    clampVariables      = clampVariables,
    markovChainParams   = markovChainParams)

def createPythonROCcfg(signalFile, bkgFiles, outFolder, csvOutFolder, treeName,
                       branchName, labels, legendPosition):
  return jinja2.Template(pythonROCcfgTemplate).render(
    signalFile     = signalFile,
    bkgFiles       = bkgFiles,
    outFolder      = outFolder,
    csvOutFolder   = csvOutFolder,
    treeName       = treeName,
    branchName     = branchName,
    labels         = labels,
    legendPosition = legendPosition)

def createBashCfg(inFileNameLocal, outFileNameLocal, inFileNameScratch,
                  outFileNameScratch, execName, pythonCfg, cmsswSrcDir):
  return jinja2.Template(jobTemplate).render(
    inFileNameLocal    = inFileNameLocal,
    outFileNameLocal   = outFileNameLocal,
    inFileNameScratch  = inFileNameScratch,
    outFileNameScratch = outFileNameScratch,
    execName           = execName,
    pythonCfg          = pythonCfg,
    cmsswSrcDir        = cmsswSrcDir)

def createSbatch(bashScript, logFile, outLocalFiles):
  return jinja2.Template(sbatchTemplate).render(
    zippedScriptLog = zip(bashScript, logFile, outLocalFiles))

def createMakefile(waitingScript, outFileNameLocalArray, scratchDir,
                   rocOutFileNames, rocCmd, rocCfg, inputBkgFiles, inputSignalFile):
  return jinja2.Template(makefileTemplate).render(
    waitingScript         = waitingScript,
    outFileNameLocalArray = outFileNameLocalArray,
    scratchDir            = scratchDir,
    rocOutFileNames       = rocOutFileNames,
    rocCmd                = rocCmd,
    rocCfg                = rocCfg,
    inputBkgFiles         = inputBkgFiles,
    inputSignalFile       = inputSignalFile)
