import FWCore.ParameterSet.Config as cms, os

evalMEMG = cms.PSet()
baseDir  = os.path.join(os.getenv('CMSSW_BASE'), 'src/tthAnalysis/tthMEM/data')

evalMEMG.tthMEM = cms.PSet(
  treeName         = cms.string('tree'),
  fileSet          = cms.VPSet(
    cms.PSet(
      identifier   = cms.string('tth'),
      fileName     = cms.string(os.path.join(baseDir, '2016/gen/out_3l_1tau_ttHJetToNonbb_M125_Tight_OS_central.root')),
      rleSelection = cms.string(os.path.expanduser('')),
    ),
    cms.PSet(
      identifier   = cms.string('ttz'),
      fileName     = cms.string(os.path.join(baseDir, '2016/gen/out_3l_1tau_TTZToLLNuNu_Tight_OS_central.root')),
      rleSelection = cms.string(os.path.expanduser('')),
    )
  ),
  madgraphFilename = cms.string(os.path.join(baseDir, 'param_card.dat')),
  higgsWidth       = cms.vdouble([-1.] + [10 * n for n in range(201)]),
  forceTauPairMass = cms.vdouble([-1., 125., 91.18]),
  logLevel         = cms.string('trace'),
  logPrecision     = cms.uint32(10),
  outputFileName   = cms.string(os.path.expanduser('~/sandbox/evalMEMG/out.root')),
  dumpToText       = cms.bool(False)
)
