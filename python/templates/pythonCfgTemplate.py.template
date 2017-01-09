import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
  fileNames   = cms.vstring('{{ inFileName }}'),
  maxEvents   = cms.int32({{ maxEvents }}),
  outputEvery = cms.uint32(50)
)

process.fwliteOutput = cms.PSet(
  fileName = cms.string('{{ outFileName }}')
)

process.logging   = cms.PSet(
  logLevel        = cms.string("info"),
  enableTimeStamp = cms.bool(True),
  enableLogging   = cms.bool(True)
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

