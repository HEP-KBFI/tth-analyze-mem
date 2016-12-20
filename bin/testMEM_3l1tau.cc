#include "tthAnalysis/tthMEM/interface/MEMInterface_3l1tau.h" // MEMInterface_3l1tau, MEMOutput_3l1tau
#include "tthAnalysis/tthMEM/interface/Logger.h" // Logger::

#include <iostream> // std::cout
#include <cstdlib> // EXIT_SUCCESS

#define VERY_QUICK 1
#define QUICK      1
#define MEDIUM     1

#define DISABLE_LOGGING 1


/**
 * @file
 *
 * This piece of code demonstrates, how to use the interface to MEM for 3l1tau channel.
 * @see test/test_MEMInterface_3l1tau for more tests that demonstrate different usage scenarios.
 */

int
main()
{
//- 1) --------------------------------------------------------------------------------------------
//--- Let's initialize the MEM
//--- This object must be persistent throughout the program's runtime
//--- (otherwise you'd lose a lot of computing time by initializing the MEM object)
#if DISABLE_LOGGING
  Logger::enableLogging(false); // disable internal MEM output (recommended)
#else
  Logger::setLogLevel(Logger::LogLevel::kInfo); // enable only for debugging purposes
#endif

  MEMInterface_3l1tau mem;

#if VERY_QUICK // 12s
  mem.maxCallsStartingPos = 100;
  mem.maxObjFunctionCalls = 1000;
#elif QUICK // 25s
  mem.maxCallsStartingPos = 5000;
  mem.maxObjFunctionCalls = 2000;
#elif MEDIUM // 2min
  mem.maxCallsStartingPos = 5000;
  mem.maxObjFunctionCalls = 10000;
#endif // original: 24min

  mem.initialize();
//--- If you want to modify some MEM parameters, do it before issuing initialize() member function,
//--- as we did in the above example. In actual analysis please do not modify the data members
//--- (in this example we meddle w/ these parameters so that testing wouldn't take lots of time).
//---
//--- It is strongly advised to disable logging as (w/ Logger::enableLogging()) well since
//---   a) printing to stdout is time-consuming and thus slows down the actual computation
//---   b) the output has little to no use unless you want to debug this executable
//---
//--- If you want to see progress of this computation, set the following bit to 1

//- 2) --------------------------------------------------------------------------------------------
//--- Let's construct the objects representing reconstructed objects

//--- Jet constructor takes: pT, eta, phi, mass
//--- Note that
//---  a) the jet ordering doesn't matter
//---  b) there must be at least 2 jets but no more than 3 jets
//--- NB! It is recommended to sort the jets by their CSV scores first, as this ordering
//--- has the highest probability to match w/ generator level b-quarks (this ordering is better
//--- than, say, the usual pT ordering).
  const MeasuredJet firstJet (184.900, -0.022,  0.460, 14.140);
  const MeasuredJet secondJet(286.300,  1.050, -2.045, 23.510);
  const MeasuredJet thirdJet ( 42.340, -1.028,  0.086,  7.765);
  const std::vector<MeasuredJet> jets{ firstJet, secondJet, thirdJet };

//--- Lepton constructor takes: pT, eta, phi, mass and charge (int)
  const MeasuredLepton leadingLepton   (117.000,  1.140,  2.136, 0.106,  1);
  const MeasuredLepton subLeadingLepton(107.800, -0.225, -0.157, 0.106, -1);
  const MeasuredLepton thirdLepton     ( 89.290,  1.639, -1.706, 0.106,  1);

//--- Hadronic tau constructor takes: pt, eta, phi, mass, charge (int), decayMode (int)
  const MeasuredHadronicTau tau(25.600, 0.410, 2.469, 0.855, -1, 1);

//--- Measured MET constructor takes: MET pT, MET phi, and (XX), (XY/YX) and (YY) components of
//---                                 MET covariance matrix
  const MeasuredMET met(164.900, 2.406, 671.900, 356.200, 775.100);

//- 3) --------------------------------------------------------------------------------------------
//--- Let's evaluate the MEM score for this event and print the results
  const MEMOutput_3l1tau result = mem(jets, leadingLepton, subLeadingLepton, thirdLepton, tau, met);
  std::cout << result << '\n';
//--- Expected output corresponding to different settings:
/**
 * VERY_QUICK
P(tth,h->tautau)              = 5.8570e-55 (1.7339e-55)
P(tth,h->WW)                  = -1.0000e+00 (-1.0000e+00)
P(ttz,z->tautau)              = 1.7738e-51 (2.2831e-52)
P(signal hypothesis) weighted = 5.8570e-55 (1.7339e-55)
P(bkg hypothesis) weighted    = 1.7738e-51 (2.2831e-52)
Likelihood ratio              = 0.000330 (0.128623)

* QUICK
P(tth,h->tautau)              = 8.4016e-53 (1.6246e-53)
P(tth,h->WW)                  = -1.0000e+00 (-1.0000e+00)
P(ttz,z->tautau)              = 6.6082e-51 (1.0005e-51)
P(signal hypothesis) weighted = 8.4016e-53 (1.6246e-53)
P(bkg hypothesis) weighted    = 6.6082e-51 (1.0005e-51)
Likelihood ratio              = 0.012554 (0.147621)

* MEDIUM
P(tth,h->tautau)              = 1.6883e-51 (2.6524e-52)
P(tth,h->WW)                  = -1.0000e+00 (-1.0000e+00)
P(ttz,z->tautau)              = 6.6437e-51 (1.3602e-51)
P(signal hypothesis) weighted = 1.6883e-51 (2.6524e-52)
P(bkg hypothesis) weighted    = 6.6437e-51 (1.3602e-51)
Likelihood ratio              = 0.202628 (0.130335)

* ORIGINAL
P(tth,h->tautau)              = 2.9619e-51 (5.4005e-52)
P(tth,h->WW)                  = -1.0000e+00 (-1.0000e+00)
P(ttz,z->tautau)              = 8.9517e-51 (9.6377e-52)
P(signal hypothesis) weighted = 2.9619e-51 (5.4005e-52)
P(bkg hypothesis) weighted    = 8.9517e-51 (9.6377e-52)
Likelihood ratio              = 0.248612 (0.061820)
*/

//--- Note that it is expected to use
//      mem(...)
//--- in a loop over the events
//--- Individual values can be obtained by inspecting member variables in MEMOutput_3l1tau struct
//--- The P(tth,h->WW) scores are currently set to -1, b/c it's not implemented, yet.
  return EXIT_SUCCESS;
}

