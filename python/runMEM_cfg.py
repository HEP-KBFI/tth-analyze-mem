import FWCore.ParameterSet.Config as cms

process = cms.PSet()
process.logging = cms.PSet(
  logLevel = cms.string("info"),
  enableTimeStamp = cms.bool(True),
  enableLogging = cms.bool(True)
)
