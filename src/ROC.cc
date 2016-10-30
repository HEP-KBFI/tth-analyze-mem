#include "tthAnalysis/tthMEM/interface/ROC.h"
#include "tthAnalysis/tthMEM/interface/Logger.h" // LOGERR, LOGINFO
#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h" // pow2()

#include <cstdlib> // EXIT_FAILURE, EXIT_SUCCESS
#include <memory> // std::shared_ptr<>
#include <functional> // std::less<>
#include <algorithm> // std::generate()

#include <boost/filesystem/path.hpp> // boost::filesystem::path
#include <boost/filesystem/operations.hpp> // boost::filesystem::exists(), ...
  // ... boost::filesystem::create_directories()

#include <TFile.h> // TFile
#include <TTree.h> // TTree
#include <TBranch.h> // TBranch
#include <TCanvas.h> // TCanvas
#include <TGraph.h> // TGraph
#include <TAxis.h> // SetLimits(), SetTitle()
#include <TMultiGraph.h> // TMultiGraph
#include <TString.h> // Form()
#include <TLegend.h> // TLegend

using namespace tthMEM;

ROC::ROC(const std::string & signalFileName,
         const std::vector<std::string> & bkgFileNames,
         const std::string & outFolderName,
         const std::string & treeName,
         const std::string & branchName,
         const std::vector<std::string> & labels)
  : signalFileName_(signalFileName)
  , bkgFileNames_(bkgFileNames)
  , outFolderName_(outFolderName)
  , treeName_(treeName)
  , branchName_(branchName)
  , labels_(labels)
{
  checkIfReady();
  readInput(signalFileName_, signal);
  for(unsigned i = 0; i < bkgFileNames_.size(); ++i)
  {
    std::vector<double> tmp;
    readInput(bkgFileNames_[i], tmp);
    bkgs.push_back(tmp);
  }

  wp = std::vector<double>(nofWPs);
  std::generate(
    wp.begin(), wp.end(),
    [this]
    {
      static int x = 0;
      return 1. * (x++) / (nofWPs - 1);
    }
  );

  LOGTRC << "Getting " << labels_[0] << " efficiencies";
  getROCcoords(signal, x);
  for(unsigned i = 0; i < bkgs.size(); ++i)
  {
    LOGTRC << "Getting " << labels_[i + 1] << " efficiencies";
    std::vector<double> y;
    getROCcoords(bkgs[i], y);
    ys.push_back(y);
  }

//--- default legend position
  legPos_ = std::vector<double>({ 0.15, 0.8, 0.3, 0.9 });
}

int
ROC::setLegendPosition(const std::vector<double> & legPos)
{
  if(legPos.size() != 4)
  {
    LOGERR << "Incorrect legend position dimensions: " << legPos.size() << ", "
           << "but must be 4";
    return EXIT_FAILURE;
  }
  if(legPos[0] >= legPos[2])
  {
    LOGERR << "x-coordinate of the upper right corner less than or equal to "
           << "the lower left corner: xl = " << legPos[0] << "; xu = " << legPos[2];
    return EXIT_FAILURE;
  }
  if(legPos[1] >= legPos[3])
  {
    LOGERR << "y-coordinate of the upper right corner less than or equal to "
           << "the lower left corner: yl = " << legPos[1] << "; yu = " << legPos[3];
    return EXIT_FAILURE;
  }
  legPos_ = legPos;
  return EXIT_SUCCESS;
}

bool
ROC::checkIfReady() const
{
  if(treeName_ == "")
  {
    LOGERR << "Tree name empty";
    return false;
  }
  if(branchName_ == "")
  {
    LOGERR << "Branch name empty";
    return false;
  }
  if(outFolderName_ == "")
  {
    LOGERR << "Output file name empty";
    return false;
  }
  else if(! boost::filesystem::exists(outFolderName_))
  {
    LOGWARN << "Folder '" << outFolderName_ << "' does not exist -> creating one";
    boost::filesystem::create_directories(outFolderName_);
  }

  std::vector<std::string> fileNames = bkgFileNames_;
  fileNames.push_back(signalFileName_);
  for(const std::string & fileName: fileNames)
    if(! boost::filesystem::exists(fileName))
    {
      LOGERR << "File '" << fileName << "' does not exists. Aborting.";
      return false;
    }
  return true;
}

void
ROC::readInput(const std::string & fileName,
               std::vector<double> & values)
{
  std::shared_ptr<TFile> f(new TFile(fileName.c_str(), "read"));
  std::shared_ptr<TTree> t(static_cast<TTree *>(f -> Get(treeName_.c_str())));

  double value;
  t -> SetBranchAddress(branchName_.c_str(), &value);
  const unsigned nofEvents = t -> GetEntries();
  for(unsigned i = 0; i < nofEvents; ++i)
  {
    t -> GetEntry(i);
    values.push_back(value);
  }
  t.reset();
  f -> Close();
  std::sort(values.begin(), values.end(), std::less<double>()); // ascending order
}

void
ROC::getROCcoords(const std::vector<double> & values,
                  std::vector<double> & coords)
{
  const double nofEvents = static_cast<double>(values.size());
//--- Loop over the working points (equally spaced in 0..1) and
//--- decide how many values are smaller than the working point.
  for(unsigned i = 0; i < wp.size(); ++i)
  {
    unsigned j = 0;
//--- Loop over values that are sorted in ascending order.
//--- If the n-th value is greater than or equal to the working point,
//--- stop there. Then there are n / N values smaller than the wp,
//--- which means that there are 1 - n / N values greater than the wp,
//--- (N being the total number of values considered), which in turn is
//--- equal to the efficiency (i.e. the percentage of values out of all
//--- values passing the working point).
    for(; j < values.size(); ++j)
      if(values[j] >= wp[i]) break;
    const double eff = 1. - j / nofEvents;
    LOGTRC << "Efficiency corresponding to WP = " << wp[i] << " is: " << eff;
    coords.push_back(eff);
  }
  coords.push_back(0.);
}

int
ROC::plotROC()
{
  if(! checkIfReady()) return EXIT_FAILURE;

//--- std::shared_ptr<> doesn't play well with TCanvas *, TGraph *, and TMultiGraph * :((
  typedef boost::filesystem::path path;
  std::vector<TGraph *> grs;
  for(unsigned i = 0; i < ys.size(); ++i)
  {
    const std::vector<double> & y = ys[i];
    if(x.size() != y.size())
    {
      LOGERR << "Something's off: the number of x and y coordinates of "
             << "the ROC curve between " << labels_[0] << " and " << labels_[i] << ' '
             << "are not equal";
      return EXIT_FAILURE;
    }
    TCanvas * c = new TCanvas(Form("ROC_%u", i), "MEM ROC", 200, 10, 700, 500);
    c -> SetGrid();
    c -> SetFillColor(0);

    TGraph * g = new TGraph(x.size(), x.data(), y.data());
    g -> SetLineColor(2 + i);
    g -> SetFillColor(0);
    g -> SetMarkerStyle(7);
    g -> SetTitle("MEM ROC curve");
    g -> Draw("ALP");
    g -> GetXaxis() -> SetLimits(0., 1.);
    g -> GetXaxis() -> SetTitle(Form("%s efficiency", labels_[0].c_str()));
    g -> GetYaxis() -> SetTitle(Form("%s efficiency", labels_[i + 1].c_str()));
    g -> SetMinimum(0.);
    g -> SetMaximum(1.);

    const path outFileName(
      path(outFolderName_) /
      path(std::string(Form("roc_%s_%s.pdf", labels_[0].c_str(), labels_[i + 1].c_str())))
    );
    c -> SaveAs(outFileName.c_str());
    grs.push_back(g);
  }

//--- if there is only one background, there is no need to create a combined plot of
//--- the ROC curves
  if(grs.size() > 1)
  {
    TCanvas * c = new TCanvas("ROC", "MEM ROC", 200, 10, 700, 500);
    c -> SetGrid();
    c -> SetFillColor(0);
    TMultiGraph * mg = new TMultiGraph();
    for(unsigned i = 0; i < grs.size(); ++i)
      mg -> Add(&*grs[i]);
    mg -> Draw("ALP");
    mg -> SetTitle("MEM ROC curve");
    mg -> GetXaxis() -> SetTitle(Form("signal (%s) efficiency", labels_[0].c_str()));
    mg -> GetYaxis() -> SetTitle("background efficiency");
    mg -> GetXaxis() -> SetLimits(0., 1.);
    mg -> SetMinimum(0.);
    mg -> SetMaximum(1.);

    TLegend * leg = new TLegend(legPos_[0], legPos_[1], legPos_[2], legPos_[3]);
    for(unsigned i = 0; i < grs.size(); ++i)
      leg -> AddEntry(&*grs[i], labels_[i + 1].c_str());
    leg -> SetBorderSize(0);
    leg -> Draw();

    const path outFileName(path(outFolderName_) / path("roc.pdf"));
    c -> SaveAs(outFileName.c_str());
  }

  return EXIT_SUCCESS;
}

double
ROC::getAUC(unsigned idx)
{
  if(! checkIfReady()) return -1.;
  double auc = 0.;
  const std::vector<double> & y = ys[idx];
  for(unsigned i = 0; i < x.size() - 1; ++i)
    auc += (x[i + 1] - x[i]) * (y[i] + y[i + 1]) / 2.;
//--- x-coordinates are sorted in descending order (@see getROCcoords())
  return -auc;
}

double
ROC::getOptimalCutoff(unsigned idx)
{
  if(! checkIfReady()) return -1.;
  double dmin = +1.e3;
  int currWPidx = -1.;
  const std::vector<double> & y = ys[idx];
  for(unsigned i = 0; i < x.size(); ++i)
    {
      const double d = std::sqrt(pow2(x[i] - 1) + pow2(y[i]));
      if(d < dmin)
      {
        dmin = d;
        currWPidx = i;
      }
    }
  LOGDBG << "Minimum distance = " << dmin << " at ("
         << x[currWPidx] << ';' << y[currWPidx] << ')';
  return wp[currWPidx];
}

