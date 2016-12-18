#ifndef MEMOUTPUT_3L1TAU_H
#define MEMOUTPUT_3L1TAU_H

struct MEMOutput_3l1tau
{
  MEMOutput_3l1tau();

  double prob_tth;
  double prob_tth_err;
  double prob_ttz;
  double prob_ttz_err;
  double prob_tth_h2ww;
  double prob_tth_h2ww_err;

  double sig;
  double sig_err;
  double sig_up;
  double sig_down;
  double bkg;
  double bkg_err;
  double bkg_up;
  double bkg_down;

  double lr;
  double lr_err;
  double lr_up;
  double lr_down;
};

#endif // MEMOUTPUT_3L1TAU_H
