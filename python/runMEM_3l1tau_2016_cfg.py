import FWCore.ParameterSet.Config as cms, os

useGen = True

process = cms.PSet()

baseDir = os.path.join(os.getenv("CMSSW_BASE"), "src/tthAnalysis/tthMEM/data/2016")
if useGen:
  baseDir = os.path.join(baseDir, "gen")

process.fwliteInput = cms.PSet(
  fileNames = cms.vstring(os.path.join(baseDir, "out_3l_1tau_ttHJetToNonbb_M125_Tight_OS_central.root")),
  # or use:
  # out_3l_1tau_TTZToLLNuNu_Tight_OS_central.root
  maxEvents = cms.int32(3),        # test with only three events
  outputEvery = cms.uint32(10000)  # never used
)

process.fwliteOutput = cms.PSet(
  fileName = cms.string(os.path.join("/scratch", os.getenv("USER"), "tree_MEM.root"))
)

process.logging = cms.PSet(
  logLevel = cms.string("trace"),
  enableTimeStamp = cms.bool(True),
  enableLogging = cms.bool(True)
)

process.tthMEM = cms.PSet(
  isMC                = cms.bool(True),
  is2016              = cms.bool(True),
  treeName            = cms.string("tree"),
  rleSelectionFile    = cms.string(""),  # run:lumi:evt numbers (one per line)
  pdfName             = cms.string("MSTW2008lo68cl"),
  madgraphFileName    = cms.string("tthAnalysis/tthMEM/data/param_card.dat"),
  integrationMode     = cms.string("VEGAS"),
  maxObjFunctionCalls = cms.uint32(150), # just for testing; proper figure: 100k+
  startingFromEntry   = cms.int64(0),
  debugPlots          = cms.uint32(1),   # use 0 if no debug plots needed
  forceGenLevel       = cms.bool(True),  # True = dumps gen lvl info regardless of clamping
  higgsWidth          = cms.double(-1.), # use negative number in case of default H width
  clampVariables      = cms.VPSet(
    cms.PSet( var = cms.string("bCosTheta1"),     useGen = cms.bool(False), useCfg = cms.bool(False), val = cms.double(0.0)),
    cms.PSet( var = cms.string("bPhi1"),          useGen = cms.bool(False), useCfg = cms.bool(False), val = cms.double(0.0)),
    cms.PSet( var = cms.string("bCosTheta2"),     useGen = cms.bool(False), useCfg = cms.bool(False), val = cms.double(0.0)),
    cms.PSet( var = cms.string("bPhi2"),          useGen = cms.bool(False), useCfg = cms.bool(False), val = cms.double(0.0)),
    cms.PSet( var = cms.string("z1"),             useGen = cms.bool(False), useCfg = cms.bool(False), val = cms.double(0.0)),
    cms.PSet( var = cms.string("tauPhi"),         useGen = cms.bool(False), useCfg = cms.bool(False), val = cms.double(0.0)),
    cms.PSet( var = cms.string("tauPhiInv"),      useGen = cms.bool(False), useCfg = cms.bool(False), val = cms.double(0.0)),
    cms.PSet( var = cms.string("tauMinvSquared"), useGen = cms.bool(False), useCfg = cms.bool(False), val = cms.double(0.0))
  ),
  markovChainParams   = cms.PSet(         # only read if the integrationMode is set to 'markovchain'
    mode                = cms.string("uniform"),
    nofBatches          = cms.uint32(27), # must divide (maxObjFunctionCalls * 0.9 / nofChains)
    nofChains           = cms.uint32(1),
    maxCallsStartingPos = cms.uint32(200),
    epsilon0            = cms.double(1.e-2),
    T0                  = cms.double(15.),
    nu                  = cms.double(0.71)
  ),
  analysisCuts = cms.PSet(
    jetAlgoRadius      = cms.double(25.),
    jetEta             = cms.double(2.4),
    jetPt              = cms.double(0.5),
    jetToLepton_dR     = cms.double(0.3),
    jetToLepton_relIso = cms.double(0.1),
  ),
)

"""
Legend:
var    = variable name you want to fix (and not integrate over)
useGen = use the value calculated from generator level information (can be used only if isMC is True)
useCfg = use the value provided in this configuration file
val    = the variable value (ignored if useCfg is False)

Therefore, for each entry there are four possible scenarious:
1) useGen = False, useCfg = False, isMC = any value:
    In this case the variable is sampled over by the lib specified by integrationMode
2) useGen = False, useCfg = True, isMC = any value:
    In this case the variable is fixated to the value specified by the current file,
    provided that the value is within expected limits
3) useGen = True, useCfg = False, isMC = True:
    Uses the value calculated from the generator level information
4) The entry is missing in clampVariables (effectively commented out):
    Reduces to the case 1)

Any other combination of boolean variables are not permitted and the program reading
the configuration file is expected to bail out before running the MEM.
"""
