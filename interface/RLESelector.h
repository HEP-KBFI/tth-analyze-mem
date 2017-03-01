#ifndef RLESELECTOR_H
#define RLESELECTOR_H

#include "tthAnalysis/tthMEM/interface/Traits.h" // _Traits<>, *_t
#include "tthAnalysis/tthMEM/interface/Exception.h" // throw_line_ext()
#include "tthAnalysis/tthMEM/interface/Logger.h" // LOG*

#include <boost/filesystem.hpp> // boost::filesystem::is_regular_file()
#include <boost/lexical_cast.hpp> // boost::lexical_cast<>(), boost::bad_lexical_cast
#include <boost/algorithm/string/join.hpp> // boost::algorithm::join()

#include <TBranch.h> // TBranch
#include <TTree.h> // TTree
#include <TFile.h> // TTree
#include <TString.h> // Form()

#include <unordered_map> // std::unordered_map<>
#include <fstream> // std::ifstream
#include <regex> // std::regex, std::regex_search(), std::smatch
#include <type_traits> // std::enable_if<>, std::
#include <system_error> // std::system_error

namespace std
{
  typedef vector<string> vstring;
}

namespace tthMEM
{
  /**
   * @brief Class for filtering events based on the run, lumi and event (RLE) numbers
   *
   * By default, the class won't filter anything, unless a valid file (aka whitelist)
   * containing a list of RLE numbers (one per line) of the format
   *   run:lumi:event
   * The class also has the capability to read run, lumi and event numbers from a given
   * TTree and write them to supplied TTree.
   */
  template<
    typename RunType  = UInt_t,
    typename LumiType = UInt_t,
    typename EvtType  = ULong64_t,
    typename = typename std::enable_if<
      std::is_integral<RunType> ::value &&
      std::is_integral<LumiType>::value &&
      std::is_integral<EvtType> ::value
    >::type
  >
  class RLESelector
  {
  public:
    typedef RunType  run_type;  ///< type of run number
    typedef LumiType lumi_type; ///< type of lumi number
    typedef EvtType  evt_type;  ///< type of event number

    /**
     * @brief Default constructor (sets only the branch names and
     *        initializes the regex object)
     */
    RLESelector()
      : run_key_ ("run")
      , lumi_key_("lumi")
      , evt_key_ ("evt")
      , run_ (0)
      , lumi_(0)
      , evt_ (0)
      , rle_regex("^([[:digit:]]+):([[:digit:]]+):([[:digit:]]+)$")
    {}
    /**
     * @brief Constructor that takes input TTree for reading
     * @param tree The TTree from which the RLE numbers will be read
     */
    RLESelector(TTree * tree)
      : RLESelector()
    {
      setBranches(tree);
    }
    /**
     * @brief Constructor that takes input TTree for reading and
     *        path to file name which holds whitelist of the RLE numbers
     * @param tree        The TTree from which the RLE number will be read
     * @param rleFileName The path to file holding the whitelised RLE numbers
     */
    RLESelector(TTree * tree,
                const std::string & rleFileName)
      : RLESelector(tree)
    {
      read(rleFileName);
    }
    /**
     * @brief Default destructor
     *
     * The destructor checks if there are any RLE numbers, which were whitelisted,
     * but never selected during the event selection. At low enough logging level,
     * the RLE numbers of whitelisted (but never selected) events will be printed
     * to stdout.
     */
    ~RLESelector()
    {
      if(! rle_.empty())
      {
        std::vstring unmatched;
        for(const auto & lumiMap: rle_)
          for(const auto & evtMap: lumiMap.second)
            for(const auto & entry: evtMap.second)
              if(! entry.second)
                unmatched.push_back(
                  RLESelector<>::str(lumiMap.first, evtMap.first, entry.first)
                );
        if(unmatched.size())
        {
          LOGTRC << "The following RLEs were selected but not matched:";
          for(const std::string missing: unmatched)
            LOGTRC << '\t' << missing;
        }
      }
    }

    /**
     * @brief Sets the (input) branches for reading
     * @param tree The TTree from which the RLE numbers are being read
     */
    void
    setBranches(TTree * tree)
    {
      tree -> SetBranchAddress(run_key_.c_str(),  &run_);
      tree -> SetBranchAddress(lumi_key_.c_str(), &lumi_);
      tree -> SetBranchAddress(evt_key_.c_str(),  &evt_);
    }

    /**
     * @brief Read the whitelisted RLE numbers from file
     * @param Path to file
     *
     * If the supplied argument is an empty string, it will be ignored
     * (meaning that the whitelist is empty and none of the events will
     * be filtered based on their RLE numbers). If the file name is there,
     * the whitelist will be generated, unless it throws in the following
     * occasions:
     *   1) the file does not exist
     *   2) there's a line in the file which doesn't respect RLE regex, i.e. form
     *        run:lumi:event
     *   3) it's impossible to cast a supplied RLE number to a given template type,
     *      e.g. the first number of run:lumi:event doesn't correspond to RunType
     *   4) other (unexpected) errors
     */
    void
    read(const std::string & rleFileName)
    {
      if(rleFileName.empty())
        return;
      if(! boost::filesystem::is_regular_file(rleFileName))
        throw_line_ext("RLE", TTHEXCEPTION_ERR_CODE_RLE_MISSING_FILE)
          << "No such file: " << rleFileName;
      std::ifstream rleFile;
      try
      {
        rleFile.open(rleFileName);
        for(std::string line; std::getline(rleFile, line); )
        {
          std::smatch what;
          if(std::regex_search(line, what, rle_regex))
          {
            try
            {
              const RunType  run  = boost::lexical_cast<RunType> (what[1]);
              const LumiType lumi = boost::lexical_cast<LumiType>(what[2]);
              const EvtType  evt  = boost::lexical_cast<EvtType> (what[3]);
              if(! rle_.count(run))
                rle_[run] = LumiMap();
              if(! rle_[run].count(lumi))
                rle_[run][lumi] = EvtMap();
              rle_[run][lumi][evt] = false;
            } catch(const boost::bad_lexical_cast & err)
            {
              throw_line_ext("RLE", TTHEXCEPTION_ERR_CODE_RLE_INVALID_TYPE)
                << "Could not convert line '" << line
                << "' to appropriate RLE types; reason: " << err.what();
            }
          }
          else
          {
            throw_line_ext("RLE", TTHEXCEPTION_ERR_CODE_RLE_INVALID_LINE)
              << "Failed to parse line '" << line << "' since it didn't match the RLE regex";
          }
        }
      } catch(const std::system_error & err)
      {
        throw_line_ext("RLE", TTHEXCEPTION_ERR_CODE_RLE_IO)
          << "Could not read file: " << rleFileName << "; reason: " << err.code().message();
      } catch(...)
      {
        throw_line_ext("RLE", TTHEXCEPTION_ERR_CODE_RLE_UNKNOWN)
          << "Could not read file for unknown reason: " << rleFileName;
      }
    }

    /**
     * @brief Sets up new branches for writing local member RLE numbers to a file
     * @param tree The TTree in which the RLE numbers will be stored
     */
    void
    initNewBranches(TTree * tree)
    {
      runBranch_  = tree -> Branch(
        run_key_.c_str(),  &run_,  Form("%s/%s", run_key_.c_str(),  _Traits<RunType>::TYPE_NAME)
      );
      lumiBranch_ = tree -> Branch(
        lumi_key_.c_str(), &lumi_, Form("%s/%s", lumi_key_.c_str(), _Traits<LumiType>::TYPE_NAME)
      );
      evtBranch_  = tree -> Branch(
        evt_key_.c_str(),  &evt_,  Form("%s/%s", evt_key_.c_str(),  _Traits<EvtType>::TYPE_NAME)
      );
    }

    /**
     * @brief Checks if current run, lumi and event number is in the whitelist
     * @return true,  if the whitelist is empty (meaning that all RLE numbers are whitelisted)
     *                or if current RLE number is contained in the whitelist (given that
     *                the whitelist is not empty);
     *         false, otherwise
     */
    bool
    has()
    {
      return has(run_, lumi_, evt_, false);
    }

    /**
     * @brief Checks if a given run, lumi and event number is allowed
     * @param run   The run number
     * @param lumi  The lumi number
     * @param evt   The event number
     * @param write Set local run, lumi and event member variables to the supplied valued
     * @return true,  if the whitelist is empty (meaning that all RLE numbers are whitelisted)
     *                or if a given RLE number is contained in the whitelist (given that
     *                the whitelist is not empty);
     *         false, otherwise
     */
    bool
    has(RunType  run,
        LumiType lumi,
        EvtType  evt,
        bool write = true)
    {
      if(write)
      {
        run_  = run;
        lumi_ = lumi;
        evt_  = evt;
      }

      if(rle_.empty())
        return true;

      const bool has_rle =
        rle_.count(run)            &&
        rle_[run].count(lumi)      &&
        rle_[run][lumi].count(evt)
      ;
      if(has_rle)
      {
        rle_[run][lumi][evt] = true;
        return true;
      }
      return false;
    }

    /**
     * @brief Resets run, lumi and event numbers to zero
     */
    void
    reset()
    {
      run_  = 0;
      lumi_ = 0;
      evt_  = 0;
    }

    /**
     * @brief Return current run:lumi:event number as a string
     * @return Colon-delimited run, lumi and event number
     */
    std::string
    str() const
    {
      return str(run_, lumi_, evt_);
    }

  private:
    typedef std::unordered_map<EvtType,  bool>    EvtMap;  ///< event number map
    typedef std::unordered_map<LumiType, EvtMap>  LumiMap; ///< lumi number map
    typedef std::unordered_map<RunType,  LumiMap> RunMap;  ///< run number map
    RunMap rle_; ///< map for storing whitelisted run, lumi and event numbers

    const std::string run_key_;  ///< default branch name of the run number
    const std::string lumi_key_; ///< default branch name of the lumi number
    const std::string evt_key_;  ///< default branch name of the event number

    RunType  run_;  ///< run number
    LumiType lumi_; ///< lumi number
    EvtType  evt_;  ///< event number

    TBranch * runBranch_  = nullptr; ///< (output) branch of the run number
    TBranch * lumiBranch_ = nullptr; ///< (output) branch of the lumi number
    TBranch * evtBranch_  = nullptr; ///< (output) branch of the event number

    const std::regex rle_regex; ///< regular expression for validation run:lumi:event numbers

    /**
     * @brief Converts any run, lumi and event numbers into a colon-delimited string
     * @param run  The run number
     * @param lumi The lumi number
     * @param evt  The event number
     * @return Colon-separated string, run:lumi:event
     */
    static std::string
    str(RunType  run,
        LumiType lumi,
        EvtType  evt)
    {
      const std::string run_str  = std::to_string(run);
      const std::string lumi_str = std::to_string(lumi);
      const std::string evt_str  = std::to_string(evt);
      return boost::algorithm::join(std::vstring({run_str, lumi_str, evt_str}), ":");
    }
  };
}

#endif // RLESELECTOR_H
