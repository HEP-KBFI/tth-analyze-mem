import FWCore.ParameterSet.Config as cms, os

roc = cms.PSet()

baseDir    = os.path.join(os.getenv("CMSSW_BASE"), "src/tthAnalysis/tthMEM/data")
rocBaseDir = os.path.join(baseDir, 'roc_test')

roc.tthMEM = cms.PSet(
  signalFile     = cms.string(os.path.join(rocBaseDir, 'testSignal.root')),
  bkgFiles       = cms.vstring(os.path.join(rocBaseDir, 'testBkg.root')),
  outFolder      = cms.string(os.path.join(baseDir, 'plots')),
  csvOutFolder   = cms.string(os.path.join(baseDir, 'csvs')),
  treeName       = cms.string('tree'),
  branchName     = cms.string('value'),
  labels         = cms.vstring('signal', 'background'),
  legendPosition = cms.vdouble(0.15, 0.78, 0.3, 0.88)
)
