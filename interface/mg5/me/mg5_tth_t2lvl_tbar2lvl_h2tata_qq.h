//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.4.2, 2016-06-10
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef mg5_tth_t2lvl_tbar2lvl_h2tata_qq_h
#define mg5_tth_t2lvl_tbar2lvl_h2tata_qq_h

#include "tthAnalysis/tthMEM/interface/mg5/me/mg5_tthz_t2lvl_tbar2lvl_hz2tata.h"
#include "tthAnalysis/tthMEM/interface/mg5/parameters/Parameters_sm_tth_t2lvl_tbar2lvl_h2tata.h"

//==========================================================================
// A class for calculating the matrix elements for
// Process: u u~ > t t~ h WEIGHTED<=4 @1
// *   Decay: t > b e+ ve WEIGHTED<=4
// *   Decay: t~ > b~ e- ve~ WEIGHTED<=4
// *   Decay: h > ta+ ta- WEIGHTED<=2
// Process: c c~ > t t~ h WEIGHTED<=4 @1
// *   Decay: t > b e+ ve WEIGHTED<=4
// *   Decay: t~ > b~ e- ve~ WEIGHTED<=4
// *   Decay: h > ta+ ta- WEIGHTED<=2
// Process: d d~ > t t~ h WEIGHTED<=4 @1
// *   Decay: t > b e+ ve WEIGHTED<=4
// *   Decay: t~ > b~ e- ve~ WEIGHTED<=4
// *   Decay: h > ta+ ta- WEIGHTED<=2
// Process: s s~ > t t~ h WEIGHTED<=4 @1
// *   Decay: t > b e+ ve WEIGHTED<=4
// *   Decay: t~ > b~ e- ve~ WEIGHTED<=4
// *   Decay: h > ta+ ta- WEIGHTED<=2
// Process: u u~ > t t~ h WEIGHTED<=4 @1
// *   Decay: t > b mu+ vm WEIGHTED<=4
// *   Decay: t~ > b~ e- ve~ WEIGHTED<=4
// *   Decay: h > ta+ ta- WEIGHTED<=2
// Process: c c~ > t t~ h WEIGHTED<=4 @1
// *   Decay: t > b mu+ vm WEIGHTED<=4
// *   Decay: t~ > b~ e- ve~ WEIGHTED<=4
// *   Decay: h > ta+ ta- WEIGHTED<=2
// Process: d d~ > t t~ h WEIGHTED<=4 @1
// *   Decay: t > b mu+ vm WEIGHTED<=4
// *   Decay: t~ > b~ e- ve~ WEIGHTED<=4
// *   Decay: h > ta+ ta- WEIGHTED<=2
// Process: s s~ > t t~ h WEIGHTED<=4 @1
// *   Decay: t > b mu+ vm WEIGHTED<=4
// *   Decay: t~ > b~ e- ve~ WEIGHTED<=4
// *   Decay: h > ta+ ta- WEIGHTED<=2
// Process: u u~ > t t~ h WEIGHTED<=4 @1
// *   Decay: t > b e+ ve WEIGHTED<=4
// *   Decay: t~ > b~ mu- vm~ WEIGHTED<=4
// *   Decay: h > ta+ ta- WEIGHTED<=2
// Process: c c~ > t t~ h WEIGHTED<=4 @1
// *   Decay: t > b e+ ve WEIGHTED<=4
// *   Decay: t~ > b~ mu- vm~ WEIGHTED<=4
// *   Decay: h > ta+ ta- WEIGHTED<=2
// Process: d d~ > t t~ h WEIGHTED<=4 @1
// *   Decay: t > b e+ ve WEIGHTED<=4
// *   Decay: t~ > b~ mu- vm~ WEIGHTED<=4
// *   Decay: h > ta+ ta- WEIGHTED<=2
// Process: s s~ > t t~ h WEIGHTED<=4 @1
// *   Decay: t > b e+ ve WEIGHTED<=4
// *   Decay: t~ > b~ mu- vm~ WEIGHTED<=4
// *   Decay: h > ta+ ta- WEIGHTED<=2
// Process: u u~ > t t~ h WEIGHTED<=4 @1
// *   Decay: t > b mu+ vm WEIGHTED<=4
// *   Decay: t~ > b~ mu- vm~ WEIGHTED<=4
// *   Decay: h > ta+ ta- WEIGHTED<=2
// Process: c c~ > t t~ h WEIGHTED<=4 @1
// *   Decay: t > b mu+ vm WEIGHTED<=4
// *   Decay: t~ > b~ mu- vm~ WEIGHTED<=4
// *   Decay: h > ta+ ta- WEIGHTED<=2
// Process: d d~ > t t~ h WEIGHTED<=4 @1
// *   Decay: t > b mu+ vm WEIGHTED<=4
// *   Decay: t~ > b~ mu- vm~ WEIGHTED<=4
// *   Decay: h > ta+ ta- WEIGHTED<=2
// Process: s s~ > t t~ h WEIGHTED<=4 @1
// *   Decay: t > b mu+ vm WEIGHTED<=4
// *   Decay: t~ > b~ mu- vm~ WEIGHTED<=4
// *   Decay: h > ta+ ta- WEIGHTED<=2
//--------------------------------------------------------------------------

class mg5_tth_t2lvl_tbar2lvl_h2tata_qq
  : public mg5_tthz_t2lvl_tbar2lvl_hz2tata
{
  public:

    // Constructor.
    mg5_tth_t2lvl_tbar2lvl_h2tata_qq()
    {
      for(std::size_t i = 0; i < nprocesses; ++i)
          jamp2[i] = nullptr;
    }

    // Destructor.
    virtual ~mg5_tth_t2lvl_tbar2lvl_h2tata_qq() override
    {
      for(std::size_t i = 0; i < nprocesses; ++i)
        if(jamp2[i])
        {
          delete jamp2[i];
          jamp2[i] = nullptr;
        }
    }

    // Initialize process.
    virtual void
    initProc(std::string param_card_name) override;

    // Calculate flavour-independent parts of cross section.
    virtual void
    sigmaKin() override;

    // Evaluate sigmaHat(sHat).
    virtual double
    sigmaHat() override;

    // Info on the subprocess.
    virtual std::string
    name() const override
    {
      return "[TTH] u u~ > b e+ ve b~ e- ve~ ta+ ta- (sm)";
    }

    // Set Higgs width
    virtual void
    setHiggsWidth(double higgsWidth) override;

    // Get matrix element vector
    const double *
    getMatrixElements() const override
    {
      return matrix_element;
    }

    static const int nprocesses = 2;

  private:

    // Private functions to calculate the matrix element for all subprocesses
    // Calculate wavefunctions
    virtual void
    calculate_wavefunctions(const int perm[],
                            const int hel[]) override;

    // Store the matrix element value from sigmaKin
    double matrix_element[nprocesses];

    // Color flows, used when selecting color
    double * jamp2[nprocesses];

    static const int nwavefuncs = 18;
    std::complex<double> w[nwavefuncs][18];
    static const int namplitudes = 2;
    std::complex<double> amp[namplitudes];

    double
    matrix_1_uux_ttxh_t_bepve_tx_bxemvex_h_taptam();

    // Pointer to the model parameters
    Parameters_sm_tth_t2lvl_tbar2lvl_h2tata * pars;
}; 


#endif  // mg5_tth_t2lvl_tbar2lvl_h2tata_qq_h
