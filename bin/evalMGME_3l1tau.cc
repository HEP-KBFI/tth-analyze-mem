#include <cstdlib> // EXIT_FAILURE, EXIT_SUCCESS
#include <string> // std::string, std::getline()
#include <cstring> // std::memset()
#include <vector> // std::vector<>
#include <array> // std::array<,>
#include <fstream> // std::ofstream::, std::scientific, std::ifstream
#include <iomanip> // std::setprecision()
#include <algorithm> // std::for_each()
#include <utility> // std::pair<,>
#include <sstream> // std::stringstream

#include "tthAnalysis/tthMEM/interface/Exception.h" // throw_line(), throw_line_ext(), tthMEMexception
#include "tthAnalysis/tthMEM/interface/Logger.h" // LOGERR, LOGINFO
#include "tthAnalysis/tthMEM/interface/mg5/me/mg5_tth_t2lvl_tbar2lvl_h2tata_qq.h" // mg5_tth_t2lvl_tbar2lvl_h2tata_qq
#include "tthAnalysis/tthMEM/interface/mg5/me/mg5_tth_t2lvl_tbar2lvl_h2tata.h" // mg5_tth_t2lvl_tbar2lvl_h2tata
#include "tthAnalysis/tthMEM/interface/mg5/me/mg5_ttz_t2lvl_tbar2lvl_z2tata.h" // mg5_ttz_t2lvl_tbar2lvl_z2tata
#include "tthAnalysis/tthMEM/interface/GeneratorLevelEvent_3l1tau.h" // tthMEM::GeneratorLevelEvent_3l1tau
#include "tthAnalysis/tthMEM/interface/general/constants.h" // tthMEM::constants
#include "tthAnalysis/tthMEM/interface/RLESelector.h" // RLESelector<>

#include <FWCore/ParameterSet/interface/ParameterSet.h> // edm::ParameterSet
#include <FWCore/PythonParameterSet/interface/MakeParameterSets.h> // edm::readPSetsFrom()
#include <FWCore/Utilities/interface/Exception.h> // cms::Exception

#include <boost/filesystem/operations.hpp> // boost::filesystem::
#include <boost/algorithm/string/join.hpp> // boost::algorithm::join()

#include <Rtypes.h> // UInt_t, ULong64_t
#include <TString.h> // Form()
#include <TFile.h> // TFile
#include <TTree.h> // TTree
#include <TBenchmark.h> // TBenchmark
#include <Math/VectorUtil.h> // ROOT::Math::VectorUtil::boost()

#define KEY_TTH_QQ "tth_qq"
#define KEY_TTH    "tth"
#define KEY_TTZ    "ttz"

/**
 * @file
 *
 * Evaluates MadGraph's (MG) matrix element (ME) for tth and tthz processes (3l1tau channel) by reading
 * the generator level information in 2016 MC samples. Plots the probability distributions.
 */

using namespace tthMEM;
namespace VectorUtil = ROOT::Math::VectorUtil;

struct PermutationLogger
{
  PermutationLogger()
    : worst{+1.e+10, default_str}
    , best{-1., default_str}
    , def{0., default_str}
  {}

  void
  process(double prob,
          unsigned bIdx,
          unsigned lIdx,
          unsigned nIdx)
  {
    const std::string curr_str = get_str(bIdx, lIdx, nIdx);
    if(is_default(curr_str))
      def = std::make_pair(prob, curr_str);

    if(prob < worst.first) worst = std::make_pair(prob, curr_str);
    if(prob > best.first)  best  = std::make_pair(prob, curr_str);
  }

  double
  get_default() const
  {
    return def.first;
  }

  double
  get_best() const
  {
    return best.first;
  }

  std::string
  get_best_str() const
  {
    return best.second;
  }

  static bool
  is_default(unsigned bIdx,
             unsigned lIdx,
             unsigned nIdx)
  {
    const std::string curr_str = get_str(bIdx, lIdx, nIdx);
    return is_default(curr_str);
  }

  static bool
  is_default(const std::string & curr_str)
  {
    return curr_str == default_str;
  }

  static std::string
  get_str(unsigned bIdx,
          unsigned lIdx,
          unsigned nIdx)
  {
    std::string curr_str = default_str;
    if(bIdx != 0) std::swap(curr_str[0], curr_str[3]);
    if(lIdx != 0) std::swap(curr_str[1], curr_str[4]);
    if(nIdx != 0) std::swap(curr_str[2], curr_str[5]);
    return curr_str;
  }

  friend std::ostream &
  operator<<(std::ostream & os,
             const PermutationLogger & logger)
  {
    os << "\tDEFAULT: ME = " << logger.def.first   << ", (" << logger.def.second   << ")\n"
       << "\tBEST:    ME = " << logger.best.first  << ", (" << logger.best.second  << ")\n"
       << "\tWORST:   ME = " << logger.worst.first << ", (" << logger.worst.second << ')';
    return os;
  }

private:
  static const std::string default_str;
  std::pair<double, std::string> worst;
  std::pair<double, std::string> best;
  std::pair<double, std::string> def;
};

const std::string PermutationLogger::default_str = "123456";

int
main(int argc,
     char * argv[])
{
  if(argc != 2)
    throw_line_ext(argv[0], TTHEXCEPTION_ERR_CODE_INVALID_PARAMETERSET)
      << "Usage: " << argv[0] << " <python config file>";

  const std::string evalMEMG = "evalMEMG";
  if(! edm::readPSetsFrom(argv[1]) -> existsAs<edm::ParameterSet>(evalMEMG.c_str()))
    throw_line_ext(argv[0], TTHEXCEPTION_ERR_CODE_INVALID_PARAMETERSET)
      << "No ParameterSet '" << evalMEMG << "' found in configuration file = " << argv[1];
  const edm::ParameterSet cfg = edm::readPSetsFrom(argv[1]) -> getParameter<edm::ParameterSet>(evalMEMG.c_str());

  const edm::ParameterSet cfg_tthMEM = cfg.getParameter<edm::ParameterSet>("tthMEM");
  const std::string treeName = cfg_tthMEM.getParameter<std::string>("treeName");
  const edm::VParameterSet fileSet = cfg_tthMEM.getParameter<edm::VParameterSet>("fileSet");
  const std::string madgraphFilename = cfg_tthMEM.getParameter<std::string>("madgraphFilename");
  const std::vector<double> higgsWidthVector = cfg_tthMEM.getParameter<std::vector<double>>("higgsWidth");
  const std::vector<double> forceTauPairMassVector = cfg_tthMEM.getParameter<std::vector<double>>("forceTauPairMass");
  const std::string logLevel = cfg_tthMEM.getParameter<std::string>("logLevel");
  const unsigned logPrecision = cfg_tthMEM.getParameter<unsigned>("logPrecision");
  const std::string outputFileName = cfg_tthMEM.getParameter<std::string>("outputFileName");
  const bool dumpToText = cfg_tthMEM.getParameter<bool>("dumpToText");
  const boost::filesystem::path outputFilePath(outputFileName);

  Logger::setLogLevel(logLevel);
  Logger::setFloatPrecision(logPrecision);

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

  mg5_tth_t2lvl_tbar2lvl_h2tata_qq tth_qq_me;
  mg5_tth_t2lvl_tbar2lvl_h2tata    tth_me;
  mg5_ttz_t2lvl_tbar2lvl_z2tata    ttz_me;
  tth_qq_me.initProc(madgraphFilename);
  tth_me.initProc(madgraphFilename);
  ttz_me.initProc(madgraphFilename);
  LOGINFO << "Declared the MG MEs";

  TFile * outFile = new TFile(outputFileName.c_str(), "recreate");
  const std::string prob_tth_qq_str = "prob_tth_qq";
  const std::string prob_tth_str    = "prob_tth";
  const std::string prob_ttz_str    = "prob_ttz";

  for(double higgsWidth: higgsWidthVector)
  {
    for(double forceTauPairMass: forceTauPairMassVector)
    {
      LOGINFO << "Higgs width = " << higgsWidth << "; tau pair mass = " << forceTauPairMass;
      if(higgsWidth > 0.)
      {
        tth_qq_me.setHiggsWidth(higgsWidth);
        tth_me.setHiggsWidth(higgsWidth);
        LOGINFO_S << "Setting Higgs width to " << higgsWidth;
      }

      const auto oIdx = [](unsigned idx) -> unsigned { return 1 - idx; };
      std::map<std::string, std::map<std::string, unsigned>> best_perms = {
        { KEY_TTH_QQ, { } },
        { KEY_TTH,    { } },
        { KEY_TTZ,    { } }
      };

      for(const edm::ParameterSet & pset: fileSet)
      {
        const std::string fileName = pset.getParameter<std::string>("fileName");
        const std::vector<std::string> id_vector = {
          pset.getParameter<std::string>("identifier"),
          (higgsWidth > 0. ? std::to_string(higgsWidth) : "original"),
          (forceTauPairMass > 0. ? std::to_string(forceTauPairMass) : "original")
        };
        const std::string id = boost::algorithm::join(id_vector, "_");

        RLESelector<> rle;
        const std::string rleSelectionFileName = pset.getParameter<std::string>("rleSelection");
        try
        {
          rle.read(rleSelectionFileName);
        }
        catch(const tthMEMexception & err)
        {
          LOGERR << err;
          return EXIT_FAILURE;
        }
        catch(...)
        {
          LOGERR << "Caught unrecognzed error while trying to initialize RLE selector";
          return EXIT_FAILURE;
        }

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

        RLESelector<>::run_type  run;
        RLESelector<>::lumi_type lumi;
        RLESelector<>::evt_type  event;
        tree -> SetBranchAddress("run",  &run);
        tree -> SetBranchAddress("lumi", &lumi);
        tree -> SetBranchAddress("evt",  &event);

        GeneratorLevelEvent_3l1tau evt;
        evt.setBranches(tree);
        LOGINFO << "Associated the input branches";

        outFile -> mkdir(id.c_str());
        outFile -> cd(id.c_str());
        TTree * outTree = new TTree(treeName.c_str(), id.c_str());

        double prob_tth_qq = 0.;
        double prob_tth    = 0.;
        double prob_ttz    = 0.;
        double prob_tth_qq_best = 0.;
        double prob_tth_best    = 0.;
        double prob_ttz_best    = 0.;
        outTree -> Branch(prob_tth_qq_str.c_str(), &prob_tth_qq, Form("%s/D", prob_tth_qq_str.c_str()));
        outTree -> Branch(prob_tth_str.c_str(),    &prob_tth,    Form("%s/D", prob_tth_str.c_str()));
        outTree -> Branch(prob_ttz_str.c_str(),    &prob_ttz,    Form("%s/D", prob_ttz_str.c_str()));
        outTree -> Branch(Form("%s_best", prob_tth_qq_str.c_str()), &prob_tth_qq_best, Form("%s_best/D", prob_tth_qq_str.c_str()));
        outTree -> Branch(Form("%s_best", prob_tth_str.c_str()),    &prob_tth_best,    Form("%s_best/D", prob_tth_str.c_str()));
        outTree -> Branch(Form("%s_best", prob_ttz_str.c_str()),    &prob_ttz_best,    Form("%s_best/D", prob_ttz_str.c_str()));
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
          if(! rle.has(run, lumi, event))
            continue;
          const std::string rle_str = rle.str();

          if(forceTauPairMass > 0.)
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

          const LorentzVector ttHorZ = (evt.genHorZ + evt.genTop[0] + evt.genTop[1]).p4(); // check
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

          PermutationLogger tth_qq_logger;
          PermutationLogger tth_logger;
          PermutationLogger ttz_logger;
          for(unsigned bIdx = 0; bIdx < 2; ++bIdx)
            for(unsigned lIdx = 0; lIdx < 2; ++lIdx)
              for(unsigned nIdx = 0; nIdx < 2; ++nIdx)
              {
                if(PermutationLogger::is_default(bIdx, lIdx, nIdx))
                {
                  if(evt.genLepFromTop[lIdx].charge()       != +1)
                    throw_line(argv[0]) << "First lepton sign is not correct";
                  if(evt.genLepFromTop[oIdx(lIdx)].charge() != -1)
                    throw_line(argv[0]) << "Second lepton sign is not correct";
                }

                const std::array<LorentzVector, 10> boosted
                {{
                  VectorUtil::boost(g[0],                                  boost),
                  VectorUtil::boost(g[1],                                  boost),
                  VectorUtil::boost(evt.genBQuarkFromTop[bIdx].p4(),       boost),
                  VectorUtil::boost(evt.genLepFromTop[lIdx].p4(),          boost),
                  VectorUtil::boost(evt.genNuFromTop[nIdx].p4(),           boost),
                  VectorUtil::boost(evt.genBQuarkFromTop[oIdx(bIdx)].p4(), boost),
                  VectorUtil::boost(evt.genLepFromTop[oIdx(lIdx)].p4(),    boost),
                  VectorUtil::boost(evt.genNuFromTop[oIdx(nIdx)].p4(),     boost),
                  VectorUtil::boost(evt.genTau[0].p4(),                    boost),
                  VectorUtil::boost(evt.genTau[1].p4(),                    boost)
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

                tth_qq_me.setMomenta(mgMomenta);
                tth_qq_me.sigmaKin();
                tth_me.setMomenta(mgMomenta);
                tth_me.sigmaKin();
                ttz_me.setMomenta(mgMomenta);
                ttz_me.sigmaKin();

                // there are two qq MEs
                const double curr_prob_tth_qq = tth_qq_me.getMatrixElements()[0];
                const double curr_prob_tth    = tth_me.getMatrixElements()[0];
                const double curr_prob_ttz    = ttz_me.getMatrixElements()[0];
                tth_qq_logger.process(curr_prob_tth_qq, bIdx, lIdx, nIdx);
                tth_logger.process(curr_prob_tth, bIdx, lIdx, nIdx);
                ttz_logger.process(curr_prob_ttz, bIdx, lIdx, nIdx);

                std::for_each(mgMomenta.begin(), mgMomenta.end(),
                  [](double * & d) { delete d; d = 0; });
              }
          prob_tth_qq = tth_qq_logger.get_default();
          prob_tth = tth_logger.get_default();
          prob_ttz = ttz_logger.get_default();
          prob_tth_best = tth_logger.get_best();
          prob_ttz_best = ttz_logger.get_best();
          prob_tth_qq_best = tth_qq_logger.get_best();
          outTree -> Fill();

          const std::map<std::string, std::string> best_perm_str = {
            { KEY_TTH_QQ, Form("%s/%s", id.c_str(), tth_qq_logger.get_best_str().c_str()) },
            { KEY_TTH,    Form("%s/%s", id.c_str(), tth_logger.get_best_str().c_str())    },
            { KEY_TTZ,    Form("%s/%s", id.c_str(), ttz_logger.get_best_str().c_str())    }
          };
          for(const auto & kv: best_perm_str)
            if(best_perms[kv.first].count(kv.second))
              ++best_perms[kv.first][kv.second];
            else
              best_perms[kv.first][kv.second] = 1;

          LOGDBG_S << "Evaluated tth qq ME for " << i << "th event (" << rle_str << ")\n"
                   << tth_qq_logger;
          LOGDBG_S << "Evaluated tth ME for " << i << "th event (" << rle_str << ")\n"
                   << tth_logger;
          LOGDBG_S << "Evaluated ttz ME for " << i << "th event (" << rle_str << ")\n"
                   << ttz_logger;

          if(dumpToText)
          {
            const std::vector<double> probs = {
              prob_tth_qq, prob_tth, prob_ttz
            };
            outTxtFile << rle_str << std::scientific << std::setprecision(6); // those are ,,sticky''
            for(std::size_t k = 0; k < probs.size(); ++k)
              outTxtFile << ',' << probs[k] << (k == probs.size() - 1 ? "\n" : "");
          }
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

      const auto strmtov = [](const std::map<std::string, unsigned> & m)
        -> std::vector<std::string>
      {
        std::vector<std::string> v;
        for(const auto & kv: m)
          v.push_back(Form("%s (%u)", kv.first.c_str(), kv.second));
        return v;
      };

      LOGDBG << "Best permutations: ";
      LOGDBG << "TTH qq: " << boost::algorithm::join(strmtov(best_perms[KEY_TTH_QQ]), ", ");
      LOGDBG << "TTH:    " << boost::algorithm::join(strmtov(best_perms[KEY_TTH]),    ", ");
      LOGDBG << "TTZ:    " << boost::algorithm::join(strmtov(best_perms[KEY_TTZ]),    ", ");
    }
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
