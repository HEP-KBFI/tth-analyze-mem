from collections import OrderedDict as OD

samples = OD()

samples["/ttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"] = OD([
  ("type", "mc"),
  ("sample_category", "signal"),
  ("process_name_specific", "ttHJetToNonbb_M125"),
  ("use_it", True),
  ("xsection", 0.2151)
])

samples["/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"] = OD([
  ("type", "mc"),
  ("sample_category", "TTZ"),
  ("process_name_specific", "TTZToLLNuNu"),
  ("use_it", True),
  ("xsection", 0.2529)
])
