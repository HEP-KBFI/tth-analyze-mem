from collections import OrderedDict as OD

samples = OD()

samples["/ttHToNonbb_M125_13TeV_powheg_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"] = OD([
  ("type", "mc"),
  ("sample_category", "signal"),
  ("process_name_specific", "ttHToNonbb_M125"),
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