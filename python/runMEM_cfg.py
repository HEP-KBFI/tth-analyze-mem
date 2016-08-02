import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
  fileNames = cms.vstring("/home/karl/ttHAnalysis/2016Aug01_dR03mvaTight/output_root/3l_1tau/Tight/ttHJetToNonbb_M125/out_3l_1tau_ttHJetToNonbb_M125_Tight_central.root"),
  maxEvents = cms.int32(-1),
  outputEvery = cms.uint32(10000)
)

process.fwliteOutput = cms.PSet(
  fileName = cms.string("/home/karl/testMEM/tree_MEM.root")
)

process.logging = cms.PSet(
  logLevel = cms.string("info"),
  enableTimeStamp = cms.bool(True),
  enableLogging = cms.bool(True)
)

process.tthMEM = cms.PSet(
  treeName = cms.string("tree"),
  pdfName = cms.string("MSTW2008lo68cl"),
  madgraphFileName = cms.string("tthAnalysis/tthMEM/data/param_card.dat"),
  integrationMode = cms.string("VEGAS"),
  maxObjFunctionCalls = cms.uint32(20000)
)
