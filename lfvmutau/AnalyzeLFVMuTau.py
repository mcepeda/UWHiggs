'''

Run LFV H->MuTau analysis in the mu+tau channel.

Authors: Maria Cepeda, Aaron Levine, Evan K. Friis, UW

'''

import MuTauTree
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
import glob
import os
#import FinalStateAnalysis.TagAndProbe.MuonPOGCorrections as MuonPOGCorrections
#import FinalStateAnalysis.TagAndProbe.H2TauCorrections as H2TauCorrections
import FinalStateAnalysis.TagAndProbe.PileupWeight as PileupWeight
import ROOT
import math

from math import sqrt, pi

data=bool ('true' in os.environ['isRealData'])
ZTauTau=  int(os.environ['isZTauTau'])

def deltaPhi(phi1, phi2):
  PHI = abs(phi1-phi2)
  if PHI<=pi:
      return PHI
  else:
      return 2*pi-PHI

def fullMT(met,mupt,taupt, metphi, muphi, tauphi):
	mux=mupt*math.cos(muphi)
	muy=mupt*math.sin(muphi)
        metx=met*math.cos(metphi)
        mety=met*math.sin(metphi)
        taux=taupt*math.cos(tauphi)
        tauy=taupt*math.sin(tauphi)
	full_et=met+mupt+taupt # for muon and tau I am approximating pt~et (M<<P)
	full_x=metx+mux+taux
        full_y=mety+muy+tauy
	full_mt_2 = full_et*full_et-full_x*full_x-full_y*full_y
	full_mt=0
	if (full_mt_2>0):
		full_mt= math.sqrt(full_mt_2)
	return full_mt

def fullPT(met,mupt,taupt, metphi, muphi, tauphi):
        mux=mupt*math.cos(muphi)
        muy=mupt*math.sin(muphi)
        metx=met*math.cos(metphi)
        mety=met*math.sin(metphi)
        taux=taupt*math.cos(tauphi)
        tauy=taupt*math.sin(tauphi)
        full_x=metx+mux+taux
        full_y=mety+muy+tauy
        full_pt_2 = full_x*full_x+full_y*full_y
        full_pt=0
        if (full_pt_2>0):
                full_pt= math.sqrt(full_pt_2)
        return full_pt

def collMass_type1(row):
        taupx=row.tPt*math.cos(row.tPhi)
        taupy=row.tPt*math.sin(row.tPhi)
        metpx = row.type1_pfMetEt*math.cos(row.type1_pfMetPhi)
        metpy = row.type1_pfMetEt*math.sin(row.type1_pfMetPhi)
        met = sqrt(metpx*metpx+metpy*metpy)

        METproj= abs(metpx*taupx+metpy*taupy)/row.tPt

        xth=row.tPt/(row.tPt+METproj)
        den=math.sqrt(xth)

        mass=row.m_t_Mass/den

        #print '%4.2f, %4.2f, %4.2f, %4.2f, %4.2f' %(scaleMass(row), den, xth, METproj,mass)


        return mass


pu_distributions = glob.glob(os.path.join(
#    'inputs', os.environ['jobid'], 'data_TauPlusX*pu.root'))
        'inputs', os.environ['jobidpu'], 'data_SingleMu*pu.root'))
pu_corrector = PileupWeight.PileupWeight('Asympt25ns', *pu_distributions)


def fakeWeight(row):
          if row.tDecayMode==0:
            return 0.62 # 0.81
          if row.tDecayMode==1:
             return 0.76 #0.80
          if row.tDecayMode==10:
            return 0.44 #0.44

class AnalyzeLFVMuTau(MegaBase):
    tree = 'mt/final/Ntuple'
    #tree = 'New_Tree'

    def __init__(self, tree, outfile, **kwargs):
        super(AnalyzeLFVMuTau, self).__init__(tree, outfile, **kwargs)
        # Use the cython wrapper
        self.tree = MuTauTree.MuTauTree(tree)
        self.out = outfile
        self.histograms = {}

    def begin(self):

        #names=["oldisotight","oldisoloose","oldisogg","oldisoboost","oldisovbf", "oldisoloosegg","oldisolooseboost","oldisoloosevbf",
         #     "newisotight","newisoloose","newisogg","newisoboost","newisovbf",  "newisoloosegg","newisolooseboost","newisoloosevbf",
          #    "noisogg","noisoboost","noisovbf",
          #    "noTauID","noiso"]
        names=["preselection","preselectionSS", "notIso","notIsoNotWeightedSS","notIsoSS","gg","boost","vbf","ggNotIso","boostNotIso","vbfNotIso","notIsoNotWeighted",
               "preselection0Jet", "preselection1Jet", "preselection2Jet","notIso0Jet", "notIso1Jet","notIso2Jet"]
        namesize = len(names)
	for x in range(0,namesize):


            self.book(names[x], "weight", "Event weight", 100, 0, 5)
            self.book(names[x], "GenWeight", "Gen level weight", 200000 ,-1000000, 1000000)
            self.book(names[x], "genHTT", "genHTT", 1000 ,0,1000)
 
            self.book(names[x], "rho", "Fastjet #rho", 100, 0, 25)
            self.book(names[x], "nvtx", "Number of vertices", 100, -0.5, 100.5)
            self.book(names[x], "prescale", "HLT prescale", 21, -0.5, 20.5)
   
            #self.book(names[x], "jet1Pt", "Muon  Pt", 300,0,300)
            #self.book(names[x], "jet2Pt", "Muon  Pt", 300,0,300)
            #self.book(names[x], "jet2Eta", "Muon  Pt", 200,-5,5)
            #self.book(names[x], "jet1Eta", "Muon  Pt", 200,-5,5)
            #self.book(names[x], "jet1PULoose ", "Muon  Pt", 3,-1,2)
            #self.book(names[x], "jet2PULoose ", "Muon  Pt", 3,-1,2)
            #self.book(names[x], "jet1PUMVA ", "Muon  Pt", 100,-1,1)
            #self.book(names[x], "jet2PUMVA ", "Muon  Pt", 100,-1,1)
            #self.book(names[x], "jet1PUTight ", "Muon  Pt", 3,-1,2)
            #self.book(names[x], "jet2PUTight ", "Muon  Pt", 3,-1,2)

 
            self.book(names[x], "mPt", "Muon  Pt", 300,0,300)
            self.book(names[x], "mEta", "Muon  eta", 100, -2.5, 2.5)
            self.book(names[x], "mMtToPfMet_type1", "Muon MT (PF Ty1)", 200, 0, 200)
            self.book(names[x], "mCharge", "Muon Charge", 5, -2, 2)

            self.book(names[x], "tPt", "Tau  Pt", 300,0,300)
            self.book(names[x], "tEta", "Tau  eta", 100, -2.5, 2.5)
            self.book(names[x], "tMtToPfMet_type1", "Tau MT (PF Ty1)", 200, 0, 200)
            self.book(names[x], "tCharge", "Tau  Charge", 5, -2, 2)
	    self.book(names[x], "tJetPt", "Tau Jet Pt" , 500, 0 ,500)	    
            self.book(names[x], "tMass", "Tau  Mass", 1000, 0, 10)
            self.book(names[x], "tLeadTrackPt", "Tau  LeadTrackPt", 300,0,300)

		       
            #self.book(names[x], "tAgainstElectronLoose", "tAgainstElectronLoose", 2,-0.5,1.5)
            self.book(names[x], "tAgainstElectronLooseMVA5", "tAgainstElectronLooseMVA5", 2,-0.5,1.5)
            #self.book(names[x], "tAgainstElectronMedium", "tAgainstElectronMedium", 2,-0.5,1.5)
            self.book(names[x], "tAgainstElectronMediumMVA5", "tAgainstElectronMediumMVA5", 2,-0.5,1.5)
            #self.book(names[x], "tAgainstElectronTight", "tAgainstElectronTight", 2,-0.5,1.5)
            self.book(names[x], "tAgainstElectronTightMVA5", "tAgainstElectronTightMVA5", 2,-0.5,1.5)
            self.book(names[x], "tAgainstElectronVTightMVA5", "tAgainstElectronVTightMVA5", 2,-0.5,1.5)


            #self.book(names[x], "tAgainstMuonLoose", "tAgainstMuonLoose", 2,-0.5,1.5)
            self.book(names[x], "tAgainstMuonLoose3", "tAgainstMuonLoose3", 2,-0.5,1.5)
            #self.book(names[x], "tAgainstMuonMedium", "tAgainstMuonMedium", 2,-0.5,1.5)
            #self.book(names[x], "tAgainstMuonTight", "tAgainstMuonTight", 2,-0.5,1.5)
            self.book(names[x], "tAgainstMuonTight3", "tAgainstMuonTight3", 2,-0.5,1.5)

            #self.book(names[x], "tAgainstMuonLooseMVA", "tAgainstMuonLooseMVA", 2,-0.5,1.5)
            #self.book(names[x], "tAgainstMuonMediumMVA", "tAgainstMuonMediumMVA", 2,-0.5,1.5)
            #self.book(names[x], "tAgainstMuonTightMVA", "tAgainstMuonTightMVA", 2,-0.5,1.5)

            self.book(names[x], "tDecayModeFinding", "tDecayModeFinding", 2,-0.5,1.5)
            self.book(names[x], "tDecayModeFindingNewDMs", "tDecayModeFindingNewDMs", 2,-0.5,1.5)
            self.book(names[x], "tDecayMode", "tDecayMode", 21,-0.5,20.5)


            self.book(names[x], "tByLooseCombinedIsolationDeltaBetaCorr3Hits", "tByLooseCombinedIsolationDeltaBetaCorr3Hits", 2,-0.5,1.5)
            self.book(names[x], "tByMediumCombinedIsolationDeltaBetaCorr3Hits", "tByMediumCombinedIsolationDeltaBetaCorr3Hits", 2,-0.5,1.5)
            self.book(names[x], "tByTightCombinedIsolationDeltaBetaCorr3Hits", "tByTightCombinedIsolationDeltaBetaCorr3Hits", 2,-0.5,1.5)

            self.book(names[x], "tByLooseIsolationMVA3newDMwLT", "tByLooseIsolationMVA3newDMwLT", 2,-0.5,1.5)
            self.book(names[x], "tByMediumIsolationMVA3newDMwLT", "tByMediumIsolationMVA3newDMwLT", 2,-0.5,1.5)
            self.book(names[x], "tByTightIsolationMVA3newDMwLT", "tByTightIsolationMVA3newDMwLT", 2,-0.5,1.5)
            self.book(names[x], "tByVTightIsolationMVA3newDMwLT", "tByVTightIsolationMVA3newDMwLT", 2,-0.5,1.5)
            self.book(names[x], "tByVVTightIsolationMVA3newDMwLT", "tByVVTightIsolationMVA3newDMwLT", 2,-0.5,1.5)

            self.book(names[x], "tByLooseIsolationMVA3oldDMwLT", "tByLooseIsolationMVA3oldDMwLT", 2,-0.5,1.5)
            self.book(names[x], "tByMediumIsolationMVA3oldDMwLT", "tByMediumIsolationMVA3oldDMwLT", 2,-0.5,1.5)
            self.book(names[x], "tByTightIsolationMVA3oldDMwLT", "tByTightIsolationMVA3oldDMwLT", 2,-0.5,1.5)
            self.book(names[x], "tByVTightIsolationMVA3oldDMwLT", "tByVTightIsolationMVA3oldDMwLT", 2,-0.5,1.5)
            self.book(names[x], "tByVVTightIsolationMVA3oldDMwLT", "tByVVTightIsolationMVA3oldDMwLT", 2,-0.5,1.5)

            #self.book(names[x], "tByLooseIsolationMVA3newDMwoLT", "tByLooseIsolationMVA3newDMwoLT", 2,-0.5,1.5)
            #self.book(names[x], "tByMediumIsolationMVA3newDMwoLT", "tByMediumIsolationMVA3newDMwoLT", 2,-0.5,1.5)
            #self.book(names[x], "tByTightIsolationMVA3newDMwoLT", "tByTightIsolationMVA3newDMwoLT", 2,-0.5,1.5)
            #self.book(names[x], "tByVTightIsolationMVA3newDMwoLT", "tByVTightIsolationMVA3newDMwoLT", 2,-0.5,1.5)
            #self.book(names[x], "tByVVTightIsolationMVA3newDMwoLT", "tByVVTightIsolationMVA3newDMwoLT", 2,-0.5,1.5)
   
            #self.book(names[x], "tByLooseIsolationMVA3oldDMwoLT", "tByLooseIsolationMVA3oldDMwoLT", 2,-0.5,1.5)
            #self.book(names[x], "tByMediumIsolationMVA3oldDMwoLT", "tByMediumIsolationMVA3oldDMwoLT", 2,-0.5,1.5)
            #self.book(names[x], "tByTightIsolationMVA3oldDMwoLT", "tByTightIsolationMVA3oldDMwoLT", 2,-0.5,1.5)
            #self.book(names[x], "tByVTightIsolationMVA3oldDMwoLT", "tByVTightIsolationMVA3oldDMwoLT", 2,-0.5,1.5)
            #self.book(names[x], "tByVVTightIsolationMVA3oldDMwoLT", "tByVVTightIsolationMVA3oldDMwoLT", 2,-0.5,1.5)

            self.book(names[x], 'mPixHits', 'Mu 1 pix hits', 10, -0.5, 9.5)
            self.book(names[x], 'mJetBtag', 'Mu 1 JetBtag', 100, -5.5, 9.5)
    	  
            self.book(names[x],"collMass_type1_1","collMass_type1_1",500,0,500);
            self.book(names[x],"collMass_type1_2","collMass_type1_2",500,0,500);

            self.book(names[x],"collMass_type1","collMass_type1",500,0,500);
            self.book(names[x],"fullMT_type1","fullMT_type1",500,0,500);
            self.book(names[x],"fullPT_type1","fullPT_type1",500,0,500);	    
    	    self.book(names[x], "LT", "ht", 400, 0, 400)
            self.book(names[x], "type1_pfMetEt", "Type1 MET", 200, 0, 200)
    
            self.book(names[x], "m_t_Mass", "Muon + Tau Mass", 200, 0, 200)
            self.book(names[x], "m_t_Pt", "Muon + Tau Pt", 200, 0, 200)
            self.book(names[x], "m_t_DR", "Muon + Tau DR", 100, 0, 10)
            self.book(names[x], "m_t_DPhi", "Muon + Tau DPhi", 100, 0, 4)
            self.book(names[x], "m_t_SS", "Muon + Tau SS", 5, -2, 2)
            self.book(names[x], "m_t_ToMETDPhi_Ty1", "Muon Tau DPhi to MET", 100, 0, 4)
    
            # Vetoes
            self.book(names[x], 'muVetoPt5IsoIdVtx', 'Number of extra muons', 5, -0.5, 4.5)
	    self.book(names[x], 'muVetoPt15IsoIdVtx', 'Number of extra muons', 5, -0.5, 4.5)
            self.book(names[x], 'tauVetoPt20Loose3HitsVtx', 'Number of extra taus', 5, -0.5, 4.5)
            self.book(names[x], 'eVetoMVAIso', 'Number of extra CiC tight electrons', 5, -0.5, 4.5)
   
            #self.book(names[x], 'jetVeto30PUCleanedTight', 'Number of extra jets', 5, -0.5, 4.5)
            #self.book(names[x], 'jetVeto30PUCleanedLoose', 'Number of extra jets', 5, -0.5, 4.5)
            self.book(names[x], 'jetVeto30', 'Number of extra jets', 5, -0.5, 4.5)	
	    #Isolation
	    self.book(names[x], 'mRelPFIsoDBDefault' ,'Muon Isolation', 100, 0.0,1.0)
   
 
            self.book(names[x], "mPhiMtPhi", "", 100, 0,4)
            self.book(names[x], "mPhiMETPhiType1", "", 100, 0,4)
            self.book(names[x], "tPhiMETPhiType1", "", 100, 0,4)

### vbf ###
            self.book(names[x], "vbfJetVeto30", "central jet veto for vbf", 5, -0.5, 4.5)
	    self.book(names[x], "vbfJetVeto20", "", 5, -0.5, 4.5)
	    self.book(names[x], "vbfMVA", "", 100, 0,0.5)
	    self.book(names[x], "vbfMass", "", 500,0,5000.0)
	    self.book(names[x], "vbfDeta", "", 100, -0.5,10.0)
            self.book(names[x], "vbfj1eta","",100,-2.5,2.5)
	    self.book(names[x], "vbfj2eta","",100,-2.5,2.5)
	    self.book(names[x], "vbfVispt","",100,0,200)
	    self.book(names[x], "vbfHrap","",100,0,5.0)
	    self.book(names[x], "vbfDijetrap","",100,0,5.0)
	    self.book(names[x], "vbfDphihj","",100,0,4)
            self.book(names[x], "vbfDphihjnomet","",100,0,4)
            self.book(names[x], "vbfNJets", "g", 5, -0.5, 4.5)
            #self.book(names[x], "vbfNJetsPULoose", "g", 5, -0.5, 4.5)
            #self.book(names[x], "vbfNJetsPUTight", "g", 5, -0.5, 4.5)

	     

    def fill_histos(self, row,name='gg', fakeRate=False):
        histos = self.histograms
        weight=1
        if (row.GenWeight!=1 and row.GenWeight!=0):
          weight=row.GenWeight

        if (fakeRate==True):
           weight=weight*fakeWeight(row)

        if (not data):
           weight*=pu_corrector(row.nTruePU)


        histos[name+'/weight'].Fill(weight)
        histos[name+'/GenWeight'].Fill(row.GenWeight)
        histos[name+'/genHTT'].Fill(row.genHTT)
        histos[name+'/rho'].Fill(row.rho, weight)
        histos[name+'/nvtx'].Fill(row.nvtx, weight)
        histos[name+'/prescale'].Fill(row.doubleMuPrescale, weight)
        #histos[name+'/jet1Pt'].Fill(row.jet1Pt, weight)
        #histos[name+'/jet2Pt'].Fill(row.jet2Pt, weight)
        #histos[name+'/jet2Eta'].Fill(row.jet2Eta, weight)
        #histos[name+'/jet1Eta'].Fill(row.jet1Eta, weight)
        #histos[name+'/jet1PULoose '].Fill(row.jet1PULoose , weight)
        #histos[name+'/jet2PULoose '].Fill(row.jet2PULoose , weight)
        #histos[name+'/jet1PUTight '].Fill(row.jet1PUTight , weight)
        #histos[name+'/jet2PUTight '].Fill(row.jet2PUTight , weight)
        #histos[name+'/jet1PUMVA '].Fill(row.jet1PUMVA , weight)
        #histos[name+'/jet2PUMVA '].Fill(row.jet2PUMVA , weight)

        histos[name+'/mPt'].Fill(row.mPt, weight)
        histos[name+'/mEta'].Fill(row.mEta, weight)
        histos[name+'/mMtToPfMet_type1'].Fill(row.mMtToPfMet_type1,weight)
        histos[name+'/mCharge'].Fill(row.mCharge, weight)
        histos[name+'/tPt'].Fill(row.tPt, weight)
        histos[name+'/tEta'].Fill(row.tEta, weight)
        histos[name+'/tMtToPfMet_type1'].Fill(row.tMtToPfMet_type1,weight)
        histos[name+'/tCharge'].Fill(row.tCharge, weight)
	histos[name+'/tJetPt'].Fill(row.tJetPt, weight)

        histos[name+'/tMass'].Fill(row.tMass,weight)
        histos[name+'/tLeadTrackPt'].Fill(row.tLeadTrackPt,weight)
		       
        #histos[name+'/tAgainstElectronLoose'].Fill(row.tAgainstElectronLoose,weight)
        histos[name+'/tAgainstElectronLooseMVA5'].Fill(row.tAgainstElectronLooseMVA5,weight)  
        #histos[name+'/tAgainstElectronMedium'].Fill(row.tAgainstElectronMedium,weight)     
        histos[name+'/tAgainstElectronMediumMVA5'].Fill(row.tAgainstElectronMediumMVA5,weight) 
        #histos[name+'/tAgainstElectronTight'].Fill(row.tAgainstElectronTight,weight)      
        histos[name+'/tAgainstElectronTightMVA5'].Fill(row.tAgainstElectronTightMVA5,weight)  
        histos[name+'/tAgainstElectronVTightMVA5'].Fill(row.tAgainstElectronVTightMVA5,weight) 


        #histos[name+'/tAgainstMuonLoose'].Fill(row.tAgainstMuonLoose,weight)
        histos[name+'/tAgainstMuonLoose3'].Fill(row.tAgainstMuonLoose3,weight)
        #histos[name+'/tAgainstMuonMedium'].Fill(row.tAgainstMuonMedium,weight)
        #histos[name+'/tAgainstMuonTight'].Fill(row.tAgainstMuonTight,weight)
        histos[name+'/tAgainstMuonTight3'].Fill(row.tAgainstMuonTight3,weight)

        #histos[name+'/tAgainstMuonLooseMVA'].Fill(row.tAgainstMuonLooseMVA,weight)
        #histos[name+'/tAgainstMuonMediumMVA'].Fill(row.tAgainstMuonMediumMVA,weight)
        #histos[name+'/tAgainstMuonTightMVA'].Fill(row.tAgainstMuonTightMVA,weight)

        histos[name+'/tDecayModeFinding'].Fill(row.tDecayModeFinding,weight)
        histos[name+'/tDecayModeFindingNewDMs'].Fill(row.tDecayModeFindingNewDMs,weight)
        histos[name+'/tDecayMode'].Fill(row.tDecayMode,weight)

        histos[name+'/tByLooseCombinedIsolationDeltaBetaCorr3Hits'].Fill(row.tByLooseCombinedIsolationDeltaBetaCorr3Hits,weight)
        histos[name+'/tByMediumCombinedIsolationDeltaBetaCorr3Hits'].Fill(row.tByMediumCombinedIsolationDeltaBetaCorr3Hits,weight)
        histos[name+'/tByTightCombinedIsolationDeltaBetaCorr3Hits'].Fill(row.tByTightCombinedIsolationDeltaBetaCorr3Hits,weight)

        histos[name+'/tByLooseIsolationMVA3newDMwLT'].Fill(row.tByLooseIsolationMVA3newDMwLT,weight)
        histos[name+'/tByMediumIsolationMVA3newDMwLT'].Fill(row.tByMediumIsolationMVA3newDMwLT,weight)
        histos[name+'/tByTightIsolationMVA3newDMwLT'].Fill(row.tByTightIsolationMVA3newDMwLT,weight)
        histos[name+'/tByVTightIsolationMVA3oldDMwLT'].Fill(row.tByVTightIsolationMVA3oldDMwLT,weight)
        histos[name+'/tByVVTightIsolationMVA3oldDMwLT'].Fill(row.tByVVTightIsolationMVA3oldDMwLT,weight)

        histos[name+'/tByLooseIsolationMVA3oldDMwLT'].Fill(row.tByLooseIsolationMVA3oldDMwLT,weight)
        histos[name+'/tByMediumIsolationMVA3oldDMwLT'].Fill(row.tByMediumIsolationMVA3oldDMwLT,weight)
        histos[name+'/tByTightIsolationMVA3oldDMwLT'].Fill(row.tByTightIsolationMVA3oldDMwLT,weight)
        histos[name+'/tByVTightIsolationMVA3oldDMwLT'].Fill(row.tByVTightIsolationMVA3oldDMwLT,weight)
        histos[name+'/tByVVTightIsolationMVA3oldDMwLT'].Fill(row.tByVVTightIsolationMVA3oldDMwLT,weight)

        #histos[name+'/tByLooseIsolationMVA3newDMwoLT'].Fill(row.tByLooseIsolationMVA3newDMwoLT,weight)
        #histos[name+'/tByMediumIsolationMVA3newDMwoLT'].Fill(row.tByMediumIsolationMVA3newDMwoLT,weight)
        #histos[name+'/tByTightIsolationMVA3newDMwoLT'].Fill(row.tByTightIsolationMVA3newDMwoLT,weight)
        #histos[name+'/tByVTightIsolationMVA3newDMwoLT'].Fill(row.tByVTightIsolationMVA3newDMwoLT,weight)
        #histos[name+'/tByVVTightIsolationMVA3newDMwoLT'].Fill(row.tByVVTightIsolationMVA3newDMwoLT,weight)
   
        #histos[name+'/tByLooseIsolationMVA3oldDMwoLT'].Fill(row.tByLooseIsolationMVA3oldDMwoLT,weight)
        #histos[name+'/tByMediumIsolationMVA3oldDMwoLT'].Fill(row.tByMediumIsolationMVA3oldDMwoLT,weight)
        #histos[name+'/tByTightIsolationMVA3oldDMwoLT'].Fill(row.tByTightIsolationMVA3oldDMwoLT,weight)
        #histos[name+'/tByVTightIsolationMVA3oldDMwoLT'].Fill(row.tByVTightIsolationMVA3oldDMwoLT,weight)
        #histos[name+'/tByVVTightIsolationMVA3oldDMwoLT'].Fill(row.tByVVTightIsolationMVA3oldDMwoLT,weight)

	histos[name+'/LT'].Fill(row.LT,weight)

        histos[name+'/collMass_type1_1'].Fill(row.m_t_collinearmass,weight)
        histos[name+'/collMass_type1_2'].Fill(row.t_m_collinearmass,weight)

        histos[name+'/collMass_type1'].Fill(collMass_type1(row),weight)
        histos[name+'/fullMT_type1'].Fill(fullMT(row.type1_pfMetEt,row.mPt,row.tPt,row.type1_pfMetPhi, row.mPhi, row.tPhi),weight)
        histos[name+'/fullPT_type1'].Fill(fullPT(row.type1_pfMetEt,row.mPt,row.tPt,row.type1_pfMetPhi, row.mPhi, row.tPhi),weight) 

	histos[name+'/type1_pfMetEt'].Fill(row.type1_pfMetEt,weight)

        histos[name+'/m_t_Mass'].Fill(row.m_t_Mass,weight)
        histos[name+'/m_t_Pt'].Fill(row.m_t_Pt,weight)
        histos[name+'/m_t_DR'].Fill(row.m_t_DR,weight)
        histos[name+'/m_t_DPhi'].Fill(row.m_t_DPhi,weight)
        histos[name+'/m_t_SS'].Fill(row.m_t_SS,weight)
	histos[name+'/m_t_ToMETDPhi_Ty1'].Fill(row.m_t_ToMETDPhi_Ty1,weight)

        histos[name+'/mPixHits'].Fill(row.mPixHits, weight)
        histos[name+'/mJetBtag'].Fill(row.mJetBtag, weight)

        histos[name+'/muVetoPt5IsoIdVtx'].Fill(row.muVetoPt5IsoIdVtx, weight)
        histos[name+'/muVetoPt15IsoIdVtx'].Fill(row.muVetoPt15IsoIdVtx, weight)
        histos[name+'/tauVetoPt20Loose3HitsVtx'].Fill(row.tauVetoPt20Loose3HitsVtx, weight)
        histos[name+'/eVetoMVAIso'].Fill(row.eVetoMVAIso, weight)
        histos[name+'/jetVeto30'].Fill(row.jetVeto30, weight)
        #histos[name+'/jetVeto30PUCleanedLoose'].Fill(row.jetVeto30PUCleanedLoose, weight)
        #histos[name+'/jetVeto30PUCleanedTight'].Fill(row.jetVeto30PUCleanedTight, weight)

	histos[name+'/mRelPFIsoDBDefault'].Fill(row.mRelPFIsoDBDefault, weight)
        
	histos[name+'/mPhiMtPhi'].Fill(deltaPhi(row.mPhi,row.tPhi),weight)
        histos[name+'/mPhiMETPhiType1'].Fill(deltaPhi(row.mPhi,row.type1_pfMetPhi),weight)
        histos[name+'/tPhiMETPhiType1'].Fill(deltaPhi(row.tPhi,row.type1_pfMetPhi),weight)
	histos[name+'/tDecayMode'].Fill(row.tDecayMode, weight)
	histos[name+'/vbfJetVeto30'].Fill(row.vbfJetVeto30, weight)
     	histos[name+'/vbfJetVeto20'].Fill(row.vbfJetVeto20, weight)
        histos[name+'/vbfMVA'].Fill(row.vbfMVA, weight)
        histos[name+'/vbfMass'].Fill(row.vbfMass, weight)
        histos[name+'/vbfDeta'].Fill(row.vbfDeta, weight)
        histos[name+'/vbfj1eta'].Fill(row.vbfj1eta, weight)
        histos[name+'/vbfj2eta'].Fill(row.vbfj2eta, weight)
        histos[name+'/vbfVispt'].Fill(row.vbfVispt, weight)
        histos[name+'/vbfHrap'].Fill(row.vbfHrap, weight)
        histos[name+'/vbfDijetrap'].Fill(row.vbfDijetrap, weight)
        histos[name+'/vbfDphihj'].Fill(row.vbfDphihj, weight)
        histos[name+'/vbfDphihjnomet'].Fill(row.vbfDphihjnomet, weight)
        histos[name+'/vbfNJets'].Fill(row.vbfNJets, weight)
        #histos[name+'/vbfNJetsPULoose'].Fill(row.vbfNJetsPULoose, weight)
        #histos[name+'/vbfNJetsPUTight'].Fill(row.vbfNJetsPUTight, weight)




    def presel(self, row):
        if not row.singleIsoMu20Pass:
            return False
        return True

    def kinematics(self, row):
        #if row.mPt < 30:
        #    return False
        if row.mPt < 30:
             return False
        if abs(row.mEta) >= 2.1:
            return False
        if row.tPt<30 :
            return False
        if abs(row.tEta)>=2.5 :
            return False
        return True

    def gg(self,row):
       if row.mPt < 45:    # there was a typo in the PAS!
           return False
       if deltaPhi(row.mPhi, row.tPhi) <2.7:
           return False
       if row.tPt < 35:
           return False
       if row.tMtToPfMet_type1 > 50:
           return False
       if row.jetVeto30!=0:
           return False
       return True

    def boost(self,row):
          if row.jetVeto30!=1:
            return False
          if row.mPt < 35:
                return False
          if row.tPt < 40:
                return False
          if row.tMtToPfMet_type1 > 50:
                return False
          return True

    def vbf(self,row):
        if row.tPt < 40:
                return False
        if row.tMtToPfMet_type1 > 50:
                return False
        if row.jetVeto30<2:
            return False
	if(row.vbfNJets<2):
	    return False
	if(abs(row.vbfDeta)<3.5):
	    return False
        if row.vbfMass < 550:
	    return False
        if row.vbfJetVeto30 > 0:
            return False
        return True

    def oppositesign(self,row):
	if row.mCharge*row.tCharge!=-1:
            return False
	return True

    #def obj1_id(self, row):
    #    return bool(row.mPFIDTight)  
 
    def obj1_id(self,row):
    	 return row.mIsGlobal and row.mIsPFMuon and (row.mNormTrkChi2<10) and (row.mMuonHits > 0) and (row.mMatchedStations > 1) and (row.mPVDXY < 0.02) and (row.mPVDZ < 0.5) and (row.mPixHits > 0) and (row.mTkLayersWithMeasurement > 5)

    def obj2_id(self, row):
	return  row.tAgainstElectronMediumMVA5 and row.tAgainstMuonTight3 and row.tDecayModeFinding

    def vetos(self,row):
		return  (bool (row.muVetoPt5IsoIdVtx<1) and bool (row.eVetoMVAIso<1) and bool (row.tauVetoPt20Loose3HitsVtx<1) )

    #def obj1_iso(self, row):
    #    return bool(row.mRelPFIsoDBDefault <0.12)
   
    def obj1_iso(self,row):
         return bool(row.mRelPFIsoDBDefault <0.1)

    def obj2_iso(self, row):
        return  row.tByTightCombinedIsolationDeltaBetaCorr3Hits

    def obj2_mediso(self, row):
	 return row.tByMediumCombinedIsolationDeltaBetaCorr3Hits

    def obj1_antiiso(self, row):
        return bool(row.mRelPFIsoDBDefault >0.2) 

    def obj2_looseiso(self, row):
        return row.tByLooseCombinedIsolationDeltaBetaCorr3Hits


    def obj2_newiso(self, row):
        return row.tByVVTightIsolationMVA3oldDMwoLT 

    #def obj2_newlooseiso(self, row):
    #    return  row.tByLooseIsolationMVA3oldDMwoLT


    def process(self):
        event =0
        sel=False
        for row in self.tree:
            if event!=row.evt:   # This is just to ensure we get the (Mu,Tau) with the highest Pt
                event=row.evt    # In principle the code saves all the MU+Tau posibilities, if an event has several combinations
                sel = False      # it will save them all.
            if sel==True:
                continue

            if not self.presel(row):
               continue 

            if ZTauTau==1 and not row.isZtautau:
               continue
            if ZTauTau==-1 and row.isZtautau:
               continue


            if not self.kinematics(row): 
                continue
 
            if not self.obj1_iso(row):
                continue

            if not self.obj1_id(row):
                continue

            if not self.vetos (row):
                continue

            #self.fill_histos(row,'noTauID')

            if not self.obj2_id (row):
                continue

            if not self.obj2_looseiso(row):
                continue

            #self.fill_histos(row,'noiso')

            if self.obj2_iso(row) and not self.oppositesign(row):
              self.fill_histos(row,'preselectionSS',False)

            if not self.obj2_iso(row) and not self.oppositesign(row):
              self.fill_histos(row,'notIsoSS',True)
              self.fill_histos(row,'notIsoNotWeightedSS',False)

            if self.obj2_iso(row) and self.oppositesign(row):  

              self.fill_histos(row,'preselection',False)
              if row.jetVeto30==0:
                self.fill_histos(row,'preselection0Jet',False)
              if row.jetVeto30==1:
                self.fill_histos(row,'preselection1Jet',False)
              if row.jetVeto30==2:
                self.fill_histos(row,'preselection2Jet',False)

              if self.gg(row):
                  self.fill_histos(row,'gg',False)

              if self.boost(row):
                  self.fill_histos(row,'boost',False)

              if self.vbf(row):
                  self.fill_histos(row,'vbf',False)

            if not self.obj2_iso(row) and self.oppositesign(row):
              self.fill_histos(row,'notIso',True)
              self.fill_histos(row,'notIsoNotWeighted',False)

              if row.jetVeto30==0:
                self.fill_histos(row,'notIso0Jet',True)
              if row.jetVeto30==1:
                self.fill_histos(row,'notIso1Jet',True)
              if row.jetVeto30==2:
                self.fill_histos(row,'notIso2Jet',True)

              if self.gg(row):
                  self.fill_histos(row,'ggNotIso',True)

              if self.boost(row):
                  self.fill_histos(row,'boostNotIso',True)

              if self.vbf(row):
                  self.fill_histos(row,'vbfNotIso',True)


            sel=True

    def finish(self):
        self.write_histos()
