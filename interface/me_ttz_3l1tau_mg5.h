//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.4.2, 2016-06-10
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef me_ttz_3l1tau_mg5_h
#define me_ttz_3l1tau_mg5_h

#include "tthAnalysis/tthMEM/interface/me_ttHorZ_3l1tau_mg5.h"
#include "Parameters_sm_ttz_3l1tau.h"

//==========================================================================
// A class for calculating the matrix elements for
// Process: g g > t t~ z WEIGHTED<=4 @1
// *   Decay: t > b e+ ve WEIGHTED<=4
// *   Decay: t~ > b~ e- ve~ WEIGHTED<=4
// *   Decay: z > ta+ ta- WEIGHTED<=2
// Process: g g > t t~ z WEIGHTED<=4 @1
// *   Decay: t > b mu+ vm WEIGHTED<=4
// *   Decay: t~ > b~ e- ve~ WEIGHTED<=4
// *   Decay: z > ta+ ta- WEIGHTED<=2
// Process: g g > t t~ z WEIGHTED<=4 @1
// *   Decay: t > b e+ ve WEIGHTED<=4
// *   Decay: t~ > b~ mu- vm~ WEIGHTED<=4
// *   Decay: z > ta+ ta- WEIGHTED<=2
// Process: g g > t t~ z WEIGHTED<=4 @1
// *   Decay: t > b mu+ vm WEIGHTED<=4
// *   Decay: t~ > b~ mu- vm~ WEIGHTED<=4
// *   Decay: z > ta+ ta- WEIGHTED<=2
//--------------------------------------------------------------------------

class me_ttz_3l1tau_mg5
  : public me_ttHorZ_3l1tau_mg5
{
  public:

    // Constructor.
    me_ttz_3l1tau_mg5() = default;

    // Destructor.
    virtual ~me_ttz_3l1tau_mg5() override {}

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
      return "[TTZ] g g > b e+ ve b~ e- ve~ ta+ ta- (sm)";
    }

    // Set Higgs width
    virtual void
    setHiggsWidth(double higgsWidth) override;

  private:

    // Private functions to calculate the matrix element for all subprocesses
    // Calculate wavefunctions
    virtual void
    calculate_wavefunctions(const int perm[],
                            const int hel[]) override;

    double
    matrix_1_gg_ttxz_t_bepve_tx_bxemvex_z_taptam();

    // Pointer to the model parameters
    Parameters_sm_ttz_3l1tau * pars;
}; 


#endif  // me_ttz_3l1tau_mg5_h
