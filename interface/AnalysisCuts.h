#ifndef ANALYSISCUTS_H
#define ANALYSISCUTS_H

struct AnalysisCuts
{
  AnalysisCuts()
    : jetPt(25.)
    , jetEta(2.4)
    , jetAlgoRadius(0.5)
    , jetToLepton_dR(0.3)
    , jetToLepton_relIso(0.1)
  {}

  double jetPt;
  double jetEta;
  double jetAlgoRadius;
  double jetToLepton_dR;
  double jetToLepton_relIso;
};

#endif // ANALYSISCUTS_H
