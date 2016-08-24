import logging, sys
from tthAnalysis.tthMEM.samples_3l1tau import samples
from tthAnalysis.tthMEM.createJobs import createJobs

if __name__ == '__main__':
  logging.basicConfig(
    stream = sys.stdout,
    level = logging.INFO,
    format = '%(asctime)s - %(levelname)s: %(message)s')

  createJobs(samples = samples,
             channel = "3l_1tau",
             version = "2016Aug01_dR03mvaTight",
             lepton_selections = ["Tight"],
             central_or_shifts = ["central"],
             execName = "runMEM_3l1tau",
             treeName = "tree",
             integrationMode = "VEGAS",
             maxObjFunctionCalls = 20000,
             nofIntegrationsPerJob = 200,
             lhRatioBranchName = "lhRatioNP",
             rocLegendPosition = [0.15, 0.78, 0.3, 0.88],
             debugPlots = True)


