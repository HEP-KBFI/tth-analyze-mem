#ifndef BJETTRANSFERFUNCTION_H
#define BJETTRANSFERFUNCTION_H

#include "tthAnalysis/tthMEM/interface/objects/MeasuredJet.h" // MeasuredJet
#include "tthAnalysis/tthMEM/interface/objects/MeasuredLepton.h" // MeasuredLepton
#include "tthAnalysis/tthMEM/interface/AnalysisCuts.h" // AnalysisCuts
#include "tthAnalysis/tthMEM/interface/general/auxFunctions.h" // pi()

namespace tthMEM
{
  namespace structs
  {
    /**
     * @brief Simple struct for holding light quark-related
     *        coefficients of its energy TF
     */
    struct qJetParams
    {
      qJetParams() = default;
      /**
       * @brief Simple constructor for determining the light quark-related
       *        coefficients (mean and stdev) of its energy TF
       * @param energy Energy of the light quark
       * @param eta    Pseudo-rapidity of the light quark
       */
      qJetParams(double energy,
                 double eta);

      /**
       * @brief Function for determining the light quark-related coefficients
       *        (mean and stdev) of its energy TF
       * @param energy Energy of the light quark
       * @param eta    Pseudo-rapidity of the light quark
       */
      void
      classify(double energy,
               double eta);

      double mu;    ///< Mean of the energy TF
      double sigma; ///< Standard deviation of the energy TF
    };

    /**
     * @brief Simple struct for holding b-quark related coefficients of its
     *        energy TF
     */
    struct bJetParams
    {
      bJetParams() = default;
      /**
       * @brief Simple constructor for determining the b-quark related
       *        coefficients (two means and stdevs) of its energy TF
       * @param energy Energy of the b-quark
       * @param eta    Pseudorapdity of the b-quark
       */
      bJetParams(double energy,
                 double eta);

      /**
       * @brief Function for determining the b-quark related coefficients
       *        (two means and stdevs) of its energy TF
       * @param energy Energy of the b-quark
       * @param eta    Pseudorapidity of the b-quark
       */
      void
      classify(double energy,
               double eta);

      qJetParams first, second;
      ///< structs for holding both means and standard deviations of the energy TF
    };
  }

  /**
   * @brief Values taken from AN2013_313_v7, p. 14
   */
  namespace constants
  {
    const double q_eta_m[2] = { 1.0,  1.0  }; // GeV^-1
    const double q_eta_a[2] = { 0.0,  0.0  }; // GeV
    const double q_eta_b[2] = { 1.56, 1.52 }; // GeV^0.5
    const double q_eta_c[2] = { 0.0,  0.13 }; // 1

    const double fb = 0.65;                     // 1
    const double b_eta_m [2] = {  1.0,  0.98 }; // GeV^-1 (1?)
    const double b_eta_m_[2] = {  0.94, 0.87 }; // GeV^-1 (1?)
    const double b_eta_n [2] = { -3.6, -4.3  }; // GeV
    const double b_eta_n_[2] = { -3.3, +9.1  }; // GeV
    const double b_eta_a [2] = {  5.7,  6.0  }; // GeV
    const double b_eta_a_[2] = {  6.6,  0.0  }; // GeV
    const double b_eta_b [2] = {  0.99, 1.9  }; // GeV^0.5
    const double b_eta_b_[2] = {  1.7,  1.1  }; // GeV^0.5
    const double b_eta_c [2] = {  0.0,  0.0  }; // 1
    const double b_eta_c_[2] = {  0.16, 0.23 }; // 1

    const double invSqrt2Pi = 1. / std::sqrt(2. * pi());
  }

  namespace functions
  {
    /**
     * @brief Delta function as the quark's energy TF
     * @param bEnergy     Computed b-quark energy
     * @param bEnergyReco Energy of the b-quark's reconstructed counterpart
     * @param bEta        Pseudo-rapidity of the b-quark
     * @return 1.
     *
     * It's delta function so using it has no effect. No function arguments
     * are actually used here. The arguments are still there to ensure
     * the same signature as in bJetTF and qJetTF
     */
    double
    deltaFunction(double bEnergy,
                  double bEnergyReco,
                  double bEta);

    /**
     * @brief Transfer function for the b-quark energy
     * @param bEnergy     Computed energy of the b-quark
     * @param bEnergyReco Energy of the b-quark's reconstructed counterpart
     * @param bEta        Pseudo-rapidity of the b-quark
     * @return Value of the transfer function
     *
     * It's modelled as a sum of two Gaussians (c.f. AN2013_313_v7, p. 13-14)
     */
    double
    bJetTF(double bEnergy,
           double bEnergyReco,
           double bEta);

    /**
     * @brief Transfer function for the light b-quark energy
     * @param qEnergy     Computed energy of the light quark
     * @param qEnergyReco Energy of the light quark's reconstructed counterpart
     * @param qEta        Pseudo-rapidity of the light quark
     * @return Value of the transfer function
     *
     * It's modelled as a single Gaussian (c.f. AN2013_313_v7, p. 13-14)
     */
    double
    qJetTF(double qEnergy,
           double qEnergyReco,
           double qEta);

    /**
     * @brief Gaussian PDF at given point
     * @param x     Given point where the PDF is evaluated at
     * @param mu    Mean of the PDF
     * @param sigma Standard deviation of the PDF
     * @return Value of the PDF at x
     */
    inline double
    gaussianPDF(double x,
                double mu,
                double sigma);

    /**
     * @brief Computes light quark's acceptance function
     * @param testQuark      The quark which is not reconstructed as a jet
     * @param leptons        All leptons in the system*
     * @param acceptedQuarks The quarks which are reconstructed as jets
     * @param quarkShiftIdx  Index of the accepted quark the energy of which is shifted
     * @param cuts           Assumed analysis cuts
     * @return Value of the acceptance function
     *
     * The quark acceptance function is implemented as described in AN2013_313_v7, p. 13-15.
     *
     * If quarkShiftIdx is negative, then no energy shift in the TF of an accepted quark
     * is not required; if it's positive, then this points to the accepted quark, the energy
     * of which must be shifted by the energy of the quark considered here (testQuark).
     * The index should be checked only if the value returned here is exactly 1.
     * If there are multiple quarks here which are not reconstructed as jets, all such indexes
     * should be collected before evaluating the TFs of reconstructed jets. The case where multiple
     * not reconstructed jets coincide with the same reconstructed jet should not be excluded:
     * all if multiple not recostructed quarks point to the same reconstructed quark, the energy
     * should in principle shifted by sum of such quarks.
     *
     * [*] All reconstructed leptons are treated as the true leptons; by ,,in the system''
     *     it is meant that all leptons that are considered in a given hypothesis.
     *
     * The function should be considered only for objects which are not reconstructed
     * as jets; for instance on a fully-sampled b-quark 4-momentum (at mass shell, so 3
     * integration variables). It's unclear, however, whether the function should also be applied
     * on a hadronic tau in an event which has only two reconstructed jets but the hypothesis
     * requires four jets.
     *
     * About integrating the light quark's TF: since integrating a gaussian function of the form
     * \f{equation}{
     *   a\exp\left[-\frac{(y - b)^{2}}{c^2}\right]
     * \f}
     *
     * from \f$ 0 \f$ to \f$ x \f$ is a difference between error functions,
     *
     * \f{equation}{
     *   \int_{0}^{x} a\exp\left[-\frac{(y - b)^{2}}{c^2}\right] dy =
     *   \frac{\sqrt{\pi}}{2} ac
     *     \left[
     *       \erf\left(\frac{b}{c}\right) - \erf\left(\frac{b - x}{c}\right)
     *     \right]\,,
     * \f}
     *
     * it follows that the integral of the TF of a light quark, \f$ T_{q}(E'; E, \eta) \f$, from
     * \f$ 0 \f$ to \f$ \frac{p_{T}}{\sin\theta} \f$ is
     *
     * \f{equation}{
     *   \frac{1}{2} \left\{ \erf(\alpha(E,\eta)) -
     *                       \erf(\alpha(E, \eta) -\frac{p_{T}}{\sqrt{2}\sigma_{q}(E,\eta)\sin\theta})
     *              \right\}\,,
     * \f}
     *
     * where \f$ \alpha(E, \eta) = \frac{\mu_{q}(E, \eta)}{\sqrt{2}\sigma_{q}(E, \eta)} \f$ and
     * \f$ \theta = \atan e^{-\eta} \f$. Since
     *
     * \f{equation}{
     *   \frac{1}{\sin\theta} = \left.\frac{1}{\sin\atan x}\right|_{x = e^{-\eta}} =
     *                        = \left.\frac{\sqrt{x^{2} + 1}}{x}\right|_{x = e^{-\eta}}\,,
     * \f}
     *
     * we can bypass computationally expensive trigonometric functions.
     */
    double
    qJetAcceptance(const MeasuredJet                 & testQuark,
                   const std::vector<MeasuredLepton> & leptons,
                   const std::vector<MeasuredJet>    & acceptedQuarks,
                   int                               & quarkShiftIdx,
                   const AnalysisCuts                & cuts = AnalysisCuts());
  }
}

#endif // BJETTRANSFERFUNCTION_H
