//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.4.2, 2016-06-10
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef MG5_tth_t2lvl_tbar2lvl_h2WW2lvlvlv_4_H
#define MG5_tth_t2lvl_tbar2lvl_h2WW2lvlvlv_4_H

#include "tthAnalysis/tthMEM/interface/mg5/parameters/Parameters_sm_tth_t2lvl_tbar2lvl_h2WW2lvlvlv.h"

//==========================================================================
// A class for calculating the matrix elements for
// Process: g g > t t~ h WEIGHTED<=4 @1
// *   Decay: t > w+ b WEIGHTED<=2
// *     Decay: w+ > mu+ vm WEIGHTED<=2
// *   Decay: t~ > w- b~ WEIGHTED<=2
// *     Decay: w- > mu- vm~ WEIGHTED<=2
// *   Decay: h > w+ w- WEIGHTED<=2
// *     Decay: w+ > e+ ve WEIGHTED<=2
// *     Decay: w- > e- ve~ WEIGHTED<=2
// Process: g g > t t~ h WEIGHTED<=4 @1
// *   Decay: t > w+ b WEIGHTED<=2
// *     Decay: w+ > mu+ vm WEIGHTED<=2
// *   Decay: t~ > w- b~ WEIGHTED<=2
// *     Decay: w- > e- ve~ WEIGHTED<=2
// *   Decay: h > w+ w- WEIGHTED<=2
// *     Decay: w+ > e+ ve WEIGHTED<=2
// *     Decay: w- > mu- vm~ WEIGHTED<=2
// Process: g g > t t~ h WEIGHTED<=4 @1
// *   Decay: t > w+ b WEIGHTED<=2
// *     Decay: w+ > e+ ve WEIGHTED<=2
// *   Decay: t~ > w- b~ WEIGHTED<=2
// *     Decay: w- > mu- vm~ WEIGHTED<=2
// *   Decay: h > w+ w- WEIGHTED<=2
// *     Decay: w+ > mu+ vm WEIGHTED<=2
// *     Decay: w- > e- ve~ WEIGHTED<=2
// Process: g g > t t~ h WEIGHTED<=4 @1
// *   Decay: t > w+ b WEIGHTED<=2
// *     Decay: w+ > e+ ve WEIGHTED<=2
// *   Decay: t~ > w- b~ WEIGHTED<=2
// *     Decay: w- > e- ve~ WEIGHTED<=2
// *   Decay: h > w+ w- WEIGHTED<=2
// *     Decay: w+ > mu+ vm WEIGHTED<=2
// *     Decay: w- > mu- vm~ WEIGHTED<=2
//--------------------------------------------------------------------------

class mg5_tth_t2lvl_tbar2lvl_h2WW2lvlvlv_4
{
  public:

    // Constructor.
    mg5_tth_t2lvl_tbar2lvl_h2WW2lvlvlv_4()
    {
      for(std::size_t i = 0; i < nprocesses; ++i)
        jamp2[i] = nullptr;
    }

    // Destructor.
    ~mg5_tth_t2lvl_tbar2lvl_h2WW2lvlvlv_4()
    {
      for(std::size_t i = 0; i < nprocesses; ++i)
        if(jamp2[i])
        {
          delete [] jamp2[i];
          jamp2[i] = nullptr;
        }
    }

    // Initialize process.
    virtual void
    initProc(std::string param_card_name);

    // Calculate flavour-independent parts of cross section.
    virtual void
    sigmaKin();

    // Evaluate sigmaHat(sHat).
    virtual double
    sigmaHat();

    // Info on the subprocess.
    virtual std::string
    name() const
    {
      return "g g > mu+ vm b mu- vm~ b~ e+ ve e- ve~ (sm)";
    }

    const std::vector<double> &
    getMasses() const
    {
      return mME;
    }

    // Get and set momenta for matrix element evaluation
    std::vector<double *>
    getMomenta()
    {
      return p;
    }

    void
    setMomenta(std::vector<double *> & momenta)
    {
      p = momenta;
    }

    void
    setInitial(int inid1,
               int inid2)
    {
      id1 = inid1;
      id2 = inid2;
    }

    // Get matrix element vector
    const double *
    getMatrixElements() const
    {
      return matrix_element;
    }

    // Constants for array limits
    static const int ninitial = 2; 
    static const int nexternal = 12; 
    static const int nprocesses = 1; 

  private:

    // Private functions to calculate the matrix element for all subprocesses
    // Calculate wavefunctions
    void
    calculate_wavefunctions(const int perm[],
                            const int hel[]);

    static const int nwavefuncs = 26; 
    std::complex<double> w[nwavefuncs][18]; 
    static const int namplitudes = 8; 
    std::complex<double> amp[namplitudes]; 

    double
    matrix_1_gg_ttxh_t_wpb_wp_mupvm_tx_wmbx_wm_mumvmx_h_wpwm_wp_epve_wm_emv_1();

    // Store the matrix element value from sigmaKin
    double matrix_element[nprocesses]; 

    // Color flows, used when selecting color
    double * jamp2[nprocesses]; 

    // Pointer to the model parameters
    Parameters_sm_tth_t2lvl_tbar2lvl_h2WW2lvlvlv * pars;

    // vector with external particle masses
    std::vector<double> mME;

    // vector with momenta (to be changed each event)
    std::vector<double *> p;

    // Initial particle ids
    int id1, id2; 
}; 


#endif  // MG5_tth_t2lvl_tbar2lvl_h2WW2lvlvlv_4_H
