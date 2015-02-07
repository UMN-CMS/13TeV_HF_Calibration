import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# replace 'myfile.root' with the source file you want to use
#Alex 8TeV 'file:/hdfs/cms/phedex/store/mc/Summer12_DR53X/DYToEE_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/02E2D67B-22F0-E111-B169-90E6BA442F3F.root' #used gsfElectrons
#Zach 13TeV "file:/hdfs/cms/phedex/store/relval/CMSSW_7_0_0/RelValZEE_13/GEN-SIM-DIGI-RECO/POSTLS170_V3_FastSim-v2/00000/ECC0AF5A-8498-E311-8708-02163E00EAC7.root"
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("file:/hdfs/cms/phedex/store/mc/Spring14dr/DYJetsToLL_M-50_13TeV-madgraph-pythia8-tauola_v2/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/00290ABC-1A06-E411-8A12-02163E00F173.root")
    )


# Output file # see : https://github.com/UMN-CMS/ZFinder/blob/master/ZFinder/Event/zfinder_extended_data_cfg.py: 24 - 27
process.TFileService = cms.Service("TFileService",
	fileName = cms.string("test.root") 
    )

# rho value for isolation
from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets # the 4 references the rParam = 0.4
process.kt6PFJetsForIsolation = kt4PFJets.clone(
	rParam = 0.6,
	doRhoFastjet = True,
	Rho_EtaMax = cms.double(2.5)
)

# Particle flow isolation
process.load("CommonTools.ParticleFlow.Isolation.pfElectronIsolation_cff") # process.pfElectronIsolationSequence, needed implicityl for process.stdElectronSequencePFIso
process.load("CommonTools.ParticleFlow.pfParticleSelection_cff") # process.pfParticleSelectionSequence

from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso
process.stdElectronSequencePFIso = setupPFElectronIso(process, 'gedGsfElectrons')

process.pfiso = cms.Sequence(
        process.pfParticleSelectionSequence
        + process.stdElectronSequencePFIso
        )

process.demo = cms.EDAnalyzer('calib')

process.load("RecoEgamma.EgammaHFProducers.hfEMClusteringSequence_cff") #Kevin: reselect hf electrons out of hf cluster

process.hfRecoEcalCandidate.intercept2DCut = cms.double(-99)# Kevin
process.hfRecoEcalCandidate.intercept2DSlope = cms.double(99)# Kevin
process.hfRecoEcalCandidate.e9e25Cut = cms.double(-1.0) # Jeremy - loosen the cut

process.p = cms.Path(
	process.kt6PFJetsForIsolation
	* process.hfRecoEcalCandidate #Kevin
	* process.pfiso
	
	* process.demo
)
