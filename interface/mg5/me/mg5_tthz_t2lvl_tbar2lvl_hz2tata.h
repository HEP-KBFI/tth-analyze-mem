#ifndef mg5_tthz_t2lvl_tbar2lvl_hz2tata_H
#define mg5_tthz_t2lvl_tbar2lvl_hz2tata_H

#include "tthAnalysis/tthMEM/interface/Exception.h" // throw_line_ext()

#include <vector> // std::vector<>

class mg5_tthz_t2lvl_tbar2lvl_hz2tata
{
public:
  mg5_tthz_t2lvl_tbar2lvl_hz2tata() = default;

  virtual ~mg5_tthz_t2lvl_tbar2lvl_hz2tata()
  { }

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

  // Set Higgs width
  virtual void
  setHiggsWidth(double higgsWidth) = 0;

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
  virtual const double *
  getMatrixElements() const
  {
    throw_line_ext("mg5_tthz_t2lvl_tbar2lvl_hz2tata", TTHEXCEPTION_ERR_CODE_UNDEFINED_FUNCTION)
      << "Function getMatrixElements() undefined for the base class!";
    return nullptr;
  }

  // Constants for array limits
  static const int ninitial = 2;
  static const int nexternal = 10;

protected:
  virtual void
  calculate_wavefunctions(const int perm[],
                          const int hel[]) = 0;

  // vector with external particle masses
  std::vector<double> mME;

  // vector with momenta (to be changed each event)
  std::vector<double *> p;

  // Initial particle ids
  int id1, id2;
};

#endif // mg5_tthz_t2lvl_tbar2lvl_hz2tata_H
