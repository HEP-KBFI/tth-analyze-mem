import logging, sys
from tthAnalysis.tthMEM.samples_3l1tau import samples
from tthAnalysis.tthMEM.createJobs import createJobs

if __name__ == '__main__':
  logging.basicConfig(
    stream = sys.stdout,
    level = logging.INFO,
    format = '%(asctime)s - %(levelname)s: %(message)s')

  clampVariables = [
    ("bCosTheta1",     False, False, 0.0),
    ("bPhi1",          False, False, 0.0),
    ("bCosTheta2",     False, False, 0.0),
    ("bPhi2",          False, False, 0.0),
    ("z1",             False, False, 0.0),
    ("tauPhi",         False, False, 0.0),
    ("tauPhiInv",      False, False, 0.0),
    ("tauMinvSquared", False, False, 0.0)
  ]
  
  markovChainParams = {
    "mode"                : "uniform",
    "nofBatches"          : 100,     # must divide (maxObjFunctionCalls * 0.9 / nofChains)
    "nofChains"           : 1,
    "maxCallsStartingPos" : 1000000,
    "epsilon0"            : 1.e-2,
    "T0"                  : 15.,
    "nu"                  : 0.71
  }

  createJobs(samples               = samples,
             channel               = "3l_1tau",
             version               = "2016Aug01_dR03mvaTight",
             lepton_selections     = ["Tight"],
             central_or_shifts     = ["central"],
             execName              = "runMEM_3l1tau",
             treeName              = "tree",
             integrationMode       = "markovchain",
             maxObjFunctionCalls   = 100000,                  # 25k per permutation
             nofIntegrationsPerJob = 25,
             lhRatioBranchName     = "lhRatioNP",
             rocLegendPosition     = [0.15, 0.78, 0.3, 0.88],
             debugPlots            = 80,                      # every 10th event is dumped to TH1D
             higgsWidth            = -1.,                     # use negative number in case of default H width
             clampVariables        = clampVariables,
             markovChainParams     = markovChainParams)


