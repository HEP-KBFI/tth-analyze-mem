#include <cstdlib> // EXIT_FAILURE, EXIT_SUCCESS
#include <string> // std::string
#include <cstring> // std::memset()
#include <vector> // std::vector<>
#include <array> // std::array<,>
#include <fstream> // std::ofstream::, std::scientific
#include <iomanip> // std::setprecision()

#include "tthAnalysis/tthMEM/interface/Exception.h" // throw_line()
#include "tthAnalysis/tthMEM/interface/Logger.h" // LOGERR, LOGINFO
#include "tthAnalysis/tthMEM/interface/me_tth_3l1tau_mg5.h" // me_tth_3l1tau_mg5
#include "tthAnalysis/tthMEM/interface/me_ttz_3l1tau_mg5.h" // me_ttz_3l1tau_mg5
#include "tthAnalysis/tthMEM/interface/GeneratorLevelEvent_3l1tau.h" // tthMEM::GeneratorLevelEvent_3l1tau
#include "tthAnalysis/tthMEM/interface/tthMEMconstants.h" // tthMEM::constants

#include <FWCore/ParameterSet/interface/ParameterSet.h> // edm::ParameterSet
#include <FWCore/PythonParameterSet/interface/MakeParameterSets.h> // edm::readPSetsFrom()
#include <FWCore/Utilities/interface/Exception.h> // cms::Exception

#include <boost/filesystem/operations.hpp> // boost::filesystem::

#include <TString.h> // Form()
#include <TFile.h> // TFile
#include <TTree.h> // TTree
#include <TBenchmark.h> // TBenchmark
#include <Math/VectorUtil.h> // ROOT::Math::VectorUtil::boost()

/**
 * @file
 *
 * Evaluates MadGraph's (MG) matrix element (ME) for tth and tthz processes (3l1tau channel) by reading
 * the generator level information in 2016 MC samples. Plots the probability distributions.
 */

using namespace tthMEM;
namespace VectorUtil = ROOT::Math::VectorUtil;

int
main(int argc,
     char * argv[])
{
  const std::string evalMEMG = "evalMEMG";
  if(! edm::readPSetsFrom(argv[1]) -> existsAs<edm::ParameterSet>(evalMEMG.c_str()))
    throw_line("tthMEM") << "No ParameterSet '" << evalMEMG << "' found in "
                         << "configuration file = " << argv[1];
  const edm::ParameterSet cfg = edm::readPSetsFrom(argv[1]) -> getParameter<edm::ParameterSet>(evalMEMG.c_str());

  const edm::ParameterSet cfg_tthMEM = cfg.getParameter<edm::ParameterSet>("tthMEM");
  const std::string treeName = cfg_tthMEM.getParameter<std::string>("treeName");
  const edm::VParameterSet fileSet = cfg_tthMEM.getParameter<edm::VParameterSet>("fileSet");
  const std::string madgraphFilename = cfg_tthMEM.getParameter<std::string>("madgraphFilename");
  const double higgsWidth = cfg_tthMEM.getParameter<double>("higgsWidth");
  const double forceTauPairMass = cfg_tthMEM.getParameter<double>("forceTauPairMass");
  const std::string logLevel = cfg_tthMEM.getParameter<std::string>("logLevel");
  const std::string outputFileName = cfg_tthMEM.getParameter<std::string>("outputFileName");
  const bool dumpToText = cfg_tthMEM.getParameter<bool>("dumpToText");
  const boost::filesystem::path outputFilePath(outputFileName);

  Logger::setLogLevel(logLevel);

  if(! boost::filesystem::is_regular_file(madgraphFilename))
  {
    LOGERR << "MadGraph parameter file '" << madgraphFilename << "' does not exist";
    return EXIT_FAILURE;
  }

  const boost::filesystem::path outputFileParentPath(outputFilePath.parent_path());
  if(! boost::filesystem::is_directory(outputFileParentPath))
  {
    LOGWARN << "Parent path of " << outputFileName << " does not exist";
    if(boost::filesystem::is_regular_file(outputFileParentPath))
    {
      LOGERR << "Parent path of " << outputFileName << " is a file and thus the directory cannot be created";
       return EXIT_FAILURE;
    }
    boost::filesystem::create_directories(outputFileParentPath);
    LOGINFO << "Created directory " << outputFileParentPath.string();
  }

  LOGINFO << "Starting the clock";
  TBenchmark clock;
  clock.Reset();
  clock.Start(evalMEMG.c_str());

  me_tth_3l1tau_mg5 tth_me;
  me_ttz_3l1tau_mg5 ttz_me;
  tth_me.initProc(madgraphFilename);
  ttz_me.initProc(madgraphFilename);
  if(higgsWidth > 0.)
  {
    tth_me.setHiggsWidth(higgsWidth);
    LOGINFO_S << "Setting Higgs width to " << higgsWidth;
  }
  LOGINFO << "Declared the MG MEs";

  TFile * outFile = new TFile(outputFileName.c_str(), "recreate");
  const std::string prob_tth_str = "prob_tth";
  const std::string prob_ttz_str = "prob_ttz";

  for(const edm::ParameterSet & pset: fileSet)
  {
    const std::string fileName = pset.getParameter<std::string>("fileName");
    const std::string id = pset.getParameter<std::string>("identifier");
    LOGINFO << "Processing " << id;

    TFile * const file = [&]()
    {
      if(! boost::filesystem::is_regular_file(fileName))
        return static_cast<TFile *>(0);
      return new TFile(fileName.c_str(), "read");
    }();
    if(! file)
    {
      LOGERR << "File '" << fileName << "' does not exist or isn't a root file";
      return EXIT_FAILURE;
    }
    TTree * const tree = [&]()
    {
      if(! file -> GetListOfKeys() -> Contains(treeName.c_str()))
        return static_cast<TTree *>(0);
      return static_cast<TTree *>(file -> Get(treeName.c_str()));
    }();
    if(! tree)
    {
      LOGERR << "Tree '" << treeName << "' does not exist in file '" << fileName << '\'';
      return EXIT_FAILURE;
    }
    LOGINFO << "Read in file '" << fileName << "' and using tree '" << treeName << '\'';

    GeneratorLevelEvent_3l1tau evt;
    evt.setBranches(tree);
    LOGINFO << "Associated the input branches";

    outFile -> mkdir(id.c_str());
    outFile -> cd(id.c_str());
    TTree * outTree = new TTree(treeName.c_str(), id.c_str());

    double prob_tth = 0;
    double prob_ttz = 0;
    outTree -> Branch(prob_tth_str.c_str(), &prob_tth, Form("%s/D", prob_tth_str.c_str()));
    outTree -> Branch(prob_ttz_str.c_str(), &prob_ttz, Form("%s/D", prob_ttz_str.c_str()));
    LOGINFO << "Associated the output branches";

    std::ofstream outTxtFile;
    const std::string outTxtFileName = (outputFileParentPath / boost::filesystem::path(id)).string() + ".csv";
    if(dumpToText)
      outTxtFile.open(outTxtFileName, std::ofstream::out);

    const unsigned long nofEvents = tree -> GetEntries();
    for(unsigned long i = 0; i < nofEvents; ++i)
    {
      tree -> GetEntry(i);
      evt.initialize();

      if(forceTauPairMass > 0)
      {
//--- Force the mass of the tau pair to a value specified in the configuration file.
//--- Since the old mass is
//---   m = sqrt(E^2 - |p|^2), where E = E_1 + E_2 and p = p_1 + p_2, and the new mass is
//---   m' = alpha * m = alpha * sqrt(E^2 - |p|^2) = sqrt((alpha * E)^2 - |alpha * p|^2),
//--- it follows that the 4-momentum of both tau leptons are multiplied by alpha, which is
//--- the ratio of new mass (m') over the original mass (m). Since we only change the mass,
//--- the 3-momentum itself remains the same for the H/Z particle.
        const double alpha = forceTauPairMass / evt.genHorZ.mass();
        for(std::size_t j = 0; j < 2; ++j)
          evt.genTau[j] = MeasuredLepton(alpha * evt.genTau[j].p4(), evt.genTau[j].charge());
        evt.genHorZ = evt.genTau[0] + evt.genTau[1];
      }

      const LorentzVector ttHorZ = (evt.genHorZ + evt.genTop[0] + evt.genTop[1]).p4();
      const double xa = (ttHorZ.e() + ttHorZ.pz()) * constants::invSqrtS;
      const double xb = (ttHorZ.e() - ttHorZ.pz()) * constants::invSqrtS;
      LOGDBG << "xa = " << xa << "; xb = " << xb;
      if(xa <= 0. || xa >= 1. || xb <= 0. || xb >= 1.)
      {
        LOGWARN << "xa or xb have unphysical values";
        prob_tth = -1.;
        prob_ttz = -1.;
        outTree -> Fill();
        continue;
      }
      const std::array<LorentzVector, 2> g
      {{
        {0., 0., +0.5 * xa * constants::sqrtS, 0.5 * xa * constants::sqrtS},
        {0., 0., -0.5 * xb * constants::sqrtS, 0.5 * xb * constants::sqrtS}
      }};
      const Vector boost(-ttHorZ.px() / ttHorZ.e(), -ttHorZ.py() / ttHorZ.e(), 0.);

      // 0, g
      // 1, g
      // 2, b, associated W+ => associated lepton +
      // 3, e+, associated W+ => associated lepton +
      // 4, ve, associated W+ => associated lepton +
      // 5, b~, associated W- => associated lepton -
      // 6, e-, associated W- => associated lepton -
      // 7, ve~, associated W- => associated lepton -
      // 8, ta+ => associated lepton +
      // 9, ta- => associated lepton -
      const std::array<LorentzVector, 10> boosted
      {{
        VectorUtil::boost(g[0],                         boost),
        VectorUtil::boost(g[1],                         boost),
        VectorUtil::boost(evt.genBQuarkFromTop[0].p4(), boost),
        VectorUtil::boost(evt.genLepFromTop[0].p4(),    boost),
        VectorUtil::boost(evt.genNuFromTop[0].p4(),     boost),
        VectorUtil::boost(evt.genBQuarkFromTop[1].p4(), boost),
        VectorUtil::boost(evt.genLepFromTop[1].p4(),    boost),
        VectorUtil::boost(evt.genNuFromTop[1].p4(),     boost),
        VectorUtil::boost(evt.genTau[0].p4(),           boost),
        VectorUtil::boost(evt.genTau[1].p4(),           boost)
      }};
      std::vector<double *> mgMomenta;
      for(const LorentzVector & lv: boosted)
      {
        double * d = new double[4];
        d[0] = lv.e();
        d[1] = lv.px();
        d[2] = lv.py();
        d[3] = lv.pz();
        mgMomenta.push_back(d);
      }
      tth_me.setMomenta(mgMomenta);
      tth_me.sigmaKin();
      ttz_me.setMomenta(mgMomenta);
      ttz_me.sigmaKin();
      prob_tth = tth_me.getMatrixElements()[0];
      prob_ttz = ttz_me.getMatrixElements()[0];
      outTree -> Fill();

      LOGDBG_S << "Evaluating tth ME for " << i << "th event: p = " << prob_tth;
      LOGDBG_S << "Evaluating ttz ME for " << i << "th event: p = " << prob_ttz;

      if(dumpToText)
        outTxtFile << std::scientific << std::setprecision(6)
                   << prob_tth << ',' << prob_ttz << '\n';
    }

    outTree -> Write();
    if(dumpToText)
    {
      outTxtFile.flush();
      outTxtFile.close();
      LOGINFO << "Wrote file: " << outTxtFileName;
    }

    if(file)
      delete file;
    if(outTree)
      delete outTree;
  }

  outFile -> Write();
  LOGINFO << "Wrote output file: " << outputFileName;

  clock.Stop(evalMEMG.c_str());
  LOGINFO << "Real time:\t" << clock.GetRealTime(evalMEMG.c_str()) << '\t'
          << "CPU time:\t" << clock.GetCpuTime(evalMEMG.c_str());

  if(outFile)
    delete outFile;

  return EXIT_SUCCESS;
}
