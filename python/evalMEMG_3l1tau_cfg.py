import FWCore.ParameterSet.Config as cms, os

evalMEMG = cms.PSet()
baseDir  = os.path.join(os.getenv('CMSSW_BASE'), 'src/tthAnalysis/tthMEM/data')

evalMEMG.tthMEM = cms.PSet(
  treeName         = cms.string('tree'),
  fileSet          = cms.VPSet(
    cms.PSet(
      identifier = cms.string('tth'),
      fileName   = cms.string(os.path.join(baseDir, '2016/out_3l_1tau_ttHToNonbb_M125_Tight_Tight_dR03mvaTight_1t0e0m0j_OS_central.root'))
    ),
    cms.PSet(
      identifier = cms.string('ttz'),
      fileName   = cms.string(os.path.join(baseDir, '2016/out_3l_1tau_TTZToLLNuNu_Tight_Tight_dR03mvaTight_1t0e0m0j_OS_central.root'))
    )
  ),
  madgraphFilename = cms.string(os.path.join(baseDir, 'param_card.dat')),
  higgsWidth       = cms.double(-1.),
  forceTauPairMass = cms.double(-1.),
  logLevel         = cms.string('info'),
  outputFileName   = cms.string(os.path.expanduser('~/sandbox/evalMEMG/out.root')),
  dumpToText       = cms.bool(True)
)