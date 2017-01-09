import FWCore.ParameterSet.Config as cms

roc = cms.PSet()

roc.tthMEM = cms.PSet(
  signalFile     = cms.string('{{ signalFile }}'),
  bkgFiles       = cms.vstring('{{ bkgFiles|join("\', \'") }}'),
  outFolder      = cms.string('{{ outFolder }}'),
  csvOutFolder   = cms.string('{{ csvOutFolder }}'),
  treeName       = cms.string('{{ treeName }}'),
  branchName     = cms.string('{{ branchName }}'),
  labels         = cms.vstring('{{ labels|join("\', \'") }}'),
  legendPosition = cms.vdouble({{ legendPosition|join(', ') }})
)

