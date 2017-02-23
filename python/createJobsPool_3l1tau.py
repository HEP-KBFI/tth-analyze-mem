#!/usr/bin/env python

import argparse, sys, logging, os, copy, jinja2, uuid
from tthAnalysis.tthMEM.samples_3l1tau import samples
from tthAnalysis.tthMEM.JobCreator import JobCreator

def getMEMversionNr(year, version, template, valueRange):
  version_nr = 1
  while True:
    hasMatch = True
    for value in valueRange:
      mem_version = template % (value, version_nr)
      if not JobCreator.isMEMPathAvailable(year, version, mem_version):
        hasMatch = False
    if not hasMatch:
      version_nr += 1
    else:
      break
  return version_nr

def createMasterMakefile(subDirs, makeFiles, masterMakeFile):
  templatesDir           = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'templates')
  masterMakefileTemplate = open(os.path.join(templatesDir, 'Makefile_master.template')).read()
  with open(masterMakeFile, 'w') as f:
    f.write(jinja2.Template(masterMakefileTemplate).render(makefiles = zip(subDirs, makeFiles)))

if __name__ == '__main__':
  logging.basicConfig(
    stream = sys.stdout,
    level  = logging.INFO,
    format = '%(asctime)s - %(levelname)s: %(message)s'
  )

  # set the help description width to 45 characters
  class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
      if text.startswith('R|'):
        return text[2:].splitlines()
      return argparse.HelpFormatter._split_lines(self, text, width)

  parser = argparse.ArgumentParser(formatter_class = lambda prog: SmartFormatter(prog, max_help_position = 50))
  parser.add_argument('-s', '--study-type', metavar = 'analysis', required = True, type = str, default = '',
                      choices = ('higgs-width', 'mxmc-initial-sweep', 'nof-calls', 'clamping', 'integrators'),
                      help = 'R|Type of analysis to run')
  parser.add_argument('-b', '--basedir', metavar = 'path', required = True, type = str, default = '',
                      help = 'R|Path to the basedir where the MEM setup and results will be stored')
  parser.add_argument('-i', '--integrator', metavar = 'integrator', required = False, type = str, default = 'markovchain',
                      choices = ('vegas', 'vamp', 'markovchain'))
  parser.add_argument('-v', '--verbose', dest = 'verbose', action = 'store_true', default = False,
                      help = 'R|Enable verbose printout')
  args = parser.parse_args()

  if args.verbose:
    logging.getLogger().setLevel(logging.DEBUG)

  if args.study_type == 'mxmc-initial-sweep' and args.integrator != 'markovchain':
    raise parser.error('Cannot use -s mxmc-initial-sweep option with other integrator than markovchain!')

  # construct a dictionary with default arguments
  clampVariables = [
    ("bCosTheta1",     False, False, 0.0),
    ("bPhi1",          False, False, 0.0),
    ("bCosTheta2",     False, False, 0.0),
    ("bPhi2",          False, False, 0.0),
    ("z1",             False, False, 0.0),
    ("tauPhi",         False, False, 0.0),
    ("tauPhiInv",      False, False, 0.0),
    ("tauMinvSquared", False, False, 0.0),
  ]
  markovChainParams = {
    "mode"                : "uniform",
    "nofBatches"          : 100,  # must divide (maxObjFunctionCalls * 0.9 / nofChains)
    "nofChains"           : 1,
    "maxCallsStartingPos" : 1000,
    "epsilon0"            : 1.e-2,
    "T0"                  : 15.,
    "nu"                  : 0.71,
  }
  defaultArguments = {
    'samples'               : samples,
    'channel'               : '3l_1tau',
    'year'                  : '2016',
    'version'               : os.path.basename(args.basedir),
    'memBaseDir'            : 'mem',
    'central_or_shifts'     : ["central"],
    'charge_selections'     : ["OS"],
    'lepton_selections'     : ["Tight"],
    'execName'              : 'runMEM_3l1tau',
    'treeName'              : 'tree',
    'rleSelectionFile'      : '',
    'integrationMode'       : args.integrator,
    'maxObjFunctionCalls'   : 100000,
    'nofIntegrationsPerJob' : 1,
    'lhRatioBranchName'     : 'lhRatioNP',
    'debugPlots'            : 0,
    'forceGenLevel'         : False,
    'higgsWidth'            : -1.,
    'clampVariables'        : clampVariables,
    'markovChainParams'     : markovChainParams,
    'comment'               : '',
    'priority'              : 'main',
    'limit'                 : 1000,
    'maxRetries'            : 3,
  }

  masterMakeFile = ''
  makeFiles, subDirs = [], []

  # use GUID as sbatch comment which is periodically checked by squeue
  # by each job instance in the job pool
  sbatchComment = str(uuid.uuid4())

  if args.study_type == 'higgs-width':
    # just vary higgsWidth parameter in some predefined range
    # let's get the version number first (so that all jobs belong to one version number)
    higgsWidths = [ -1., ]
    higgsWidths.extend([i / 100. for i in range(1, 10)])
    higgsWidths.extend([i / 10.  for i in range(1, 21)])
    hw_template = "mem_higgsWidth_%.3fGeV_v%d"
    version_nr = getMEMversionNr(
      defaultArguments['year'], defaultArguments['version'], hw_template, higgsWidths
    )

    # now create a bunch of jobs, one for each Higgs width
    for higgsWidth in higgsWidths:
      higgsArgs = copy.deepcopy(defaultArguments)
      higgsWidth_comment = '%.3f GeV' % higgsWidth if higgsWidth > 0. else 'original'
      higgsArgs['higgsWidth'] = higgsWidth
      higgsArgs['comment']    = "Higgs width %s" % higgsWidth_comment
      higgsArgs['memBaseDir'] = hw_template % (higgsWidth, version_nr)
      jobs = JobCreator(**higgsArgs)
      jobs.createJobs(sbatchComment)

      subDirs.append(jobs.memDir)
      makeFiles.append(jobs.makeFile)
      if not masterMakeFile:
        masterMakeFile = os.path.join(jobs.baseDir, "Makefile_higgsWidth_v%d" % version_nr)

      logging.debug("Created jobs for Higgs width of %.3f GeV (version %d)" %
                    (higgsWidth, version_nr))

  elif args.study_type == 'mxmc-initial-sweep':
    # just vary markovChainParams's variable maxCallsStartingPos
    maxCallsStartingPositions = [100000, 25000, 10000, 2500, 1000]
    mcsp_template = "mem_maxCallsStartingPos_%d_v%d"
    version_nr = getMEMversionNr(
      defaultArguments['year'], defaultArguments['version'], mcsp_template, maxCallsStartingPositions
    )

    # now create a bunch of jobs, one for each max call starting position
    for maxCallsStartingPos in maxCallsStartingPositions:
      mcspArgs = copy.deepcopy(defaultArguments)
      mcspArgs['markovChainParams']['maxCallsStartingPos'] = maxCallsStartingPos
      mcspArgs['comment']                                  = "Max calls starting pos %d" % maxCallsStartingPos
      mcspArgs['memBaseDir']                               = mcsp_template % (maxCallsStartingPos, version_nr)
      jobs = JobCreator(**mcspArgs)
      jobs.createJobs(sbatchComment)

      subDirs.append(jobs.memDir)
      makeFiles.append(jobs.makeFile)
      if not masterMakeFile:
        masterMakeFile = os.path.join(jobs.baseDir, "Makefile_MCSP_v%d" % version_nr)

      logging.debug("Created jobs for max calls starting position of %.3f (version %d)" %
                    (maxCallsStartingPos, version_nr))

  elif args.study_type == 'nof-calls':
    # vary maxObjFunctionCalls
    #                    1M    500k    200k    100k    50k    25k    10k    5k
    nofCallsList = [1000000, 500000, 200000, 100000, 50000, 25000, 10000, 5000]
    nofCalls_template = "mem_nofCalls_{integrator}_%d_v%d".format(integrator = defaultArguments['integrationMode'])
    version_nr = getMEMversionNr(
      defaultArguments['year'], defaultArguments['version'], nofCalls_template, nofCallsList
    )

    # now create a bunch of jobs, one for each max call starting position
    for nofCalls in nofCallsList:
      nofCallsArgs = copy.deepcopy(defaultArguments)
      nofCallsArgs['maxObjFunctionCalls'] = nofCalls
      nofCallsArgs['comment']             = "NOF calls %d" % nofCalls
      nofCallsArgs['memBaseDir']          = nofCalls_template % (nofCalls, version_nr)
      jobs = JobCreator(**nofCallsArgs)
      jobs.createJobs(sbatchComment)

      subDirs.append(jobs.memDir)
      makeFiles.append(jobs.makeFile)
      if not masterMakeFile:
        masterMakeFile = os.path.join(
          jobs.baseDir, "Makefile_nofCalls_{integrator}_v%d".format(
            integrator = defaultArguments['integrationMode']
          ) % version_nr
        )

      logging.debug("Created jobs for integrator %s at nof calls of %d (version %d)" %
                    (defaultArguments['integrationMode'], nofCalls, version_nr))

  elif args.study_type == 'clamping':
    # we want to test the following scenarios:
    # 1) gen lvl objects, all clamped
    # 2) gen lvl objects, all but z1 clamped
    # 3) gen lvl objects, only z1 clamped
    # 4) gen lvl objects, none clamped
    # 5) reco lvl objects, none clamped
    # therefore, the naming scheme is the following:
    #   mem_clamp_[gen|reco]_[all|none|z1|not-z1]_v$(version number)
    clampList = [ "reco_none", "gen_none", "gen_all", "gen_z1", "gen_all-but-z1" ]
    clamp_template = "mem_clamp_%s_v%d"
    version_nr = getMEMversionNr(
      defaultArguments['year'], defaultArguments['version'], clamp_template, clampList
    )

    # now create a bunch of jobs, one for each max call starting position
    for clampType in clampList:
      clampArgs = copy.deepcopy(defaultArguments)
      clampArgs['comment']    = '%s objects, %s clamped' % tuple(clampType.split('_'))
      clampArgs['memBaseDir'] = clamp_template % (clampType, version_nr)

      if clampType == "reco_none":
        pass
      elif clampType == "gen_none":
        clampArgs['forceGenLevel'] = True
      elif clampType == "gen_all":
        clampArgs['forceGenLevel'] = True
        for clampVariable_idx in range(len(clampArgs['clampVariables'])):
          clampVariable = clampArgs['clampVariables'][clampVariable_idx]
          tmpClampVariable = list(clampArgs['clampVariables'][clampVariable_idx])
          tmpClampVariable[1] = True
          clampArgs['clampVariables'][clampVariable_idx] = tuple(tmpClampVariable)
      elif clampType == "gen_z1":
        clampArgs['forceGenLevel'] = True
        for clampVariable_idx in range(len(clampArgs['clampVariables'])):
          clampVariable = clampArgs['clampVariables'][clampVariable_idx]
          if clampVariable[0] == 'z1':
            tmpClampVariable = list(clampArgs['clampVariables'][clampVariable_idx])
            tmpClampVariable[1] = True
            clampArgs['clampVariables'][clampVariable_idx] = tuple(tmpClampVariable)
      elif clampType == "gen_all-but-z1":
        clampArgs['forceGenLevel'] = True
        for clampVariable_idx in range(len(clampArgs['clampVariables'])):
          clampVariable = clampArgs['clampVariables'][clampVariable_idx]
          if clampVariable[0] != 'z1':
            tmpClampVariable = list(clampArgs['clampVariables'][clampVariable_idx])
            tmpClampVariable[1] = True
            clampArgs['clampVariables'][clampVariable_idx] = tuple(tmpClampVariable)
      else:
        raise ValueError("Internal error: no such clamping type: {clampType}".format(clampType = clampType))

      jobs = JobCreator(**clampArgs)
      jobs.createJobs(sbatchComment)

      subDirs.append(jobs.memDir)
      makeFiles.append(jobs.makeFile)
      if not masterMakeFile:
        masterMakeFile = os.path.join(jobs.baseDir, "Makefile_clamp_v%d" % version_nr)

      logging.debug("Created jobs for clamping regime '%s' (version %d)" %
                    (clampArgs['comment'], version_nr))

  elif args.study_type == 'integrators':
    # we want to test MEM with different integrators: Markov Chain integrator, VEGAS and VAMP
    # if you want to perform the test at different number of function calls, edit the defaultArguments dict manually
    integratorList = ['markovchain', 'VEGAS', 'VAMP']
    integrator_template = 'mem_integrator_%s'
    version_nr = getMEMversionNr(
      defaultArguments['year'], defaultArguments['version'], integrator_template, integratorList
    )

    # now create a bunch of jobs, one for each max call starting position
    for integratorType in integratorList:
      integratorArgs = copy.deepcopy(defaultArguments)
      integratorArgs['comment']         = '%s integrator' % integratorType
      integratorArgs['memBaseDir']      = integrator_template % integratorType
      integratorArgs['integrationMode'] = integratorType

      jobs = JobCreator(**integratorArgs)
      jobs.createJobs(sbatchComment)

      subDirs.append(jobs.memDir)
      makeFiles.append(jobs.makeFile)
      if not masterMakeFile:
        masterMakeFile = os.path.join(jobs.baseDir, 'Makefile_integrators_v%d' % version_nr)

      logging.debug("Created jobs for integrator '%s' (version %d)" %
                    (integratorType, version_nr))

  else:
    raise ValueError("Internal error: unimplemented study type: {study_type}".format(study_type = args.study_type))

  # create a ,,master'' Makefile
  assert (masterMakeFile)
  createMasterMakefile(subDirs, makeFiles, masterMakeFile)
  logging.info("Run it with:\n\tmake -f %s -j 8 &> %s.log" % (masterMakeFile, masterMakeFile))
