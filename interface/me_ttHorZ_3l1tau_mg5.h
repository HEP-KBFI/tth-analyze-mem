#ifndef ME_TTHORZ_3L1TAU_MG5_H
#define ME_TTHORZ_3L1TAU_MG5_H

#include <string> // std::string
#include <vector> // std::vector<>
#include <complex> // std::complex<>

class me_ttHorZ_3l1tau_mg5
{
public:
  me_ttHorZ_3l1tau_mg5() = default;

  virtual ~me_ttHorZ_3l1tau_mg5()
  {
    if(jamp2[0])
    {
      delete jamp2[0];
      jamp2[0] = 0;
    }
  }

  // Initialize process.
  virtual void
  initProc(std::string param_card_name) = 0;

  // Calculate flavour-independent parts of cross section.
  virtual void
  sigmaKin() = 0;

  // Evaluate sigmaHat(sHat).
  virtual double
  sigmaHat() = 0;

  // Info on the subprocess.
  virtual std::string
  name() const = 0;

  const std::vector<double> &
  getMasses() const
  {
    return mME;
  }

  // Get and set momenta for matrix element evaluation
  std::vector <double *>
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
  static const int nexternal = 10;
  static const int nprocesses = 1;

protected:
  virtual void
  calculate_wavefunctions(const int perm[],
                          const int hel[]) = 0;
  static const int nwavefuncs = 22;
  std::complex<double> w[nwavefuncs][18];
  static const int namplitudes = 8;
  std::complex<double> amp[namplitudes];

  // Store the matrix element value from sigmaKin
  double matrix_element[nprocesses];

  // Color flows, used when selecting color
  double * jamp2[nprocesses];

  // vector with external particle masses
  std::vector<double> mME;

  // vector with momenta (to be changed each event)
  std::vector <double *> p;

  // Initial particle ids
  int id1, id2;
};

#endif // ME_TTHORZ_3L1TAU_MG5_H
