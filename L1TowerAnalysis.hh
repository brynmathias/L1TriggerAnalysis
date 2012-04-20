#include "L1Ntuple.h"
#include "hist.C"
#include "Style.C"
#include "TMath.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <iostream>
class L1TowerAnalysis : public L1Ntuple
{
  public :

    //constructor
  L1TowerAnalysis(std::string filename) : L1Ntuple(filename) {}
  L1TowerAnalysis() {}
  ~L1TowerAnalysis() {}

    //main function macro : arguments can be adpated to your need
  void run(Long64_t nevents, TString outputname ,TString TriggerBit);

  private :
  bool PassTrig(int ib,int bx);
  bool MatchJet(int RecoJetIdx);
  bool PassHLT( TString TrigBit);
  bool LooseID( int Jet );
  bool L1_TrippleJet(double Thresh1, double Thresh2, double Thresh3);
  bool L1_Quad_Cen_Jet(double Threshold);
  bool L1_Double_Cen_Jet(double Threshold);
  bool L1_Double_Tau_Jet(double Threshold, double Eta);
  bool NL1Jets_Threshold(int n, double threshold);
  
  std::pair <int,int> ReturnMatchedJet(int RecoJetIdx);
  bool MatchemuJet( int Jet );
  std::pair <int,int> ReturnMatchedemuJet(int RecoJetIdx);
  double ReturnMatchedemuQuantity( std::pair<int,int> matchJet,int Quantity);
  void BookHistos();
  double ReturnMatchedQuantity( std::pair<int,int> matchJet,int Quantity);
  virtual double deltaPhi(double phi1, double phi2);
  virtual double deltaR(double eta1, double phi1, double eta2, double phi2);
  double MaxL1Et();
  int leadingOfflineJet();
//your private methods can be declared here
// histos
  TH1F *RefJets;
  TH1F *ET10;
  TH1F *l1emuTauEt;
  TH1F *L1HT;
  TH1F *L1MET;
  TH1F *AllL1Jets;
  TH1F *centralJets;
  TH1F *NonForward;
  TH2F *emu_extraJetCorrelation;
  TH1F *events;
  TH1F* forward;
  TH1F* QuadJet20;
  TH1F *RefJetsemu;
  TH1F *ET10emu;
  TH1F *l1TauEtemu;
  TH1F *l1TauEt;
  TH1F *l1emuTauEtemu;
  TH1F *L1HTemu;
  TH1F *L1METemu;
  TH1F *AllL1Jetsemu;
  TH1F *centralJetsemu;
  TH1F* forwardemu;
  TH1F* QuadJet20emu;



};

double L1TowerAnalysis::deltaPhi(double phi1, double phi2) {
  double result = phi1 - phi2;
  if(fabs(result) > 9999) return result;
  while (result > TMath::Pi()) result -= 2*TMath::Pi();
  while (result <= -TMath::Pi()) result += 2*TMath::Pi();
  return result;
}

double L1TowerAnalysis::deltaR(double eta1, double phi1, double eta2, double phi2) {
  double deta = eta1 - eta2;
  double dphi = deltaPhi(phi1, phi2);
  return sqrt(deta*deta + dphi*dphi);
}

  bool L1TowerAnalysis::NL1Jets_Threshold(int n, double threshold){
    int N = 0;
    for( unsigned int a = 0; a < l1extra_->cenJetEt.size(); a++){
      if( l1extra_->cenJetEt[a] > threshold ) N++;
    }
    for( unsigned int a = 0; a < l1extra_->tauJetEt.size(); a++){
      if( l1extra_->tauJetEt[a] > threshold ) N++;
    }
    return N >= n;
  }



bool L1TowerAnalysis::L1_Double_Cen_Jet(double Threshold){
  if (l1extra_->cenJetEt.size() + l1extra_->tauJetEt.size() < 2) return false;
  if (l1extra_->cenJetEt[0] > Threshold && l1extra_->cenJetEt[1] > Threshold) return true;
  return false;
}
bool L1TowerAnalysis::L1_Quad_Cen_Jet(double Threshold){
  if (l1extra_->cenJetEt.size() < 4) return false;
  if (l1extra_->cenJetEt[0] > Threshold && l1extra_->cenJetEt[1] > Threshold && l1extra_->cenJetEt[2] > Threshold && l1extra_->cenJetEt[3] > Threshold) return true;
  return false;
}
bool L1TowerAnalysis::L1_TrippleJet(double Thresh1, double Thresh2, double Thresh3){
  if (l1extra_->cenJetEt.size() < 4) return false;
  if (l1extra_->cenJetEt[0] > Thresh1 && l1extra_->cenJetEt[1] > Thresh2 && l1extra_->cenJetEt[2] > Thresh3) return true;
  return false;
}

bool L1TowerAnalysis::L1_Double_Tau_Jet(double Threshold, double Eta){
  if(l1extra_->tauJetEt.size() < 2 )return false;
  if((l1extra_->tauJetEt[1] > Threshold && l1extra_->tauJetEt[0] > Threshold) && ( fabs(l1extra_->tauJetEta[0]) < Eta && fabs(l1extra_->tauJetEta[1]) < Eta ))return true;
  return false;
}


bool L1TowerAnalysis::MatchJet(int RecoJetIdx){
  double minDeltaR = 999999.;
  for(unsigned int i = 0; i < l1extra_->nCenJets; i++){
    if( deltaR(recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx], l1extra_->cenJetEta[i], l1extra_->cenJetPhi[i]) < minDeltaR)
    {
        minDeltaR = deltaR(recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx], l1extra_->cenJetEta[i], l1extra_->cenJetPhi[i]);
    }
  }

  for(unsigned int i = 0; i < l1extra_->nTauJets; i++){
    if( deltaR(recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx], l1extra_->tauJetEta[i], l1extra_->tauJetPhi[i]) < minDeltaR)
    {
      minDeltaR = deltaR(recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx], l1extra_->tauJetEta[i], l1extra_->tauJetPhi[i]);
    }
  }

  for(unsigned int i = 0; i < l1extra_->nFwdJets; i++){
    if( deltaR(recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx], l1extra_->fwdJetEta[i], l1extra_->fwdJetPhi[i]) < minDeltaR)
    {
      minDeltaR = deltaR(recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx], l1extra_->fwdJetEta[i], l1extra_->fwdJetPhi[i]);
    }
  }
  if(minDeltaR < .5){ return true; }
  else return false;
}

bool L1TowerAnalysis::MatchemuJet(int RecoJetIdx){
  double minDeltaR = 999999.;
  for(unsigned int i = 0; i < l1emuextra_->nCenJets; i++){
    if( deltaR(recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx], l1emuextra_->cenJetEta[i], l1emuextra_->cenJetPhi[i]) < minDeltaR)
    {
        minDeltaR = deltaR(recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx], l1emuextra_->cenJetEta[i], l1emuextra_->cenJetPhi[i]);
    }
  }

  for(unsigned int i = 0; i < l1emuextra_->nTauJets; i++){
    if( deltaR(recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx], l1emuextra_->tauJetEta[i], l1emuextra_->tauJetPhi[i]) < minDeltaR)
    {
      minDeltaR = deltaR(recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx], l1emuextra_->tauJetEta[i], l1emuextra_->tauJetPhi[i]);
    }
  }

  for(unsigned int i = 0; i < l1emuextra_->nFwdJets; i++){
    if( deltaR(recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx], l1emuextra_->fwdJetEta[i], l1emuextra_->fwdJetPhi[i]) < minDeltaR)
    {
      minDeltaR = deltaR(recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx], l1emuextra_->fwdJetEta[i], l1emuextra_->fwdJetPhi[i]);
    }
  }
  if(minDeltaR < .5){ return true; }
  else return false;
}



std::pair<int,int> L1TowerAnalysis::ReturnMatchedJet(int RecoJetIdx){
  double minDeltaR = 9999.;
  std::pair <int,int> matchedJet(-1,-1);
  //printf("Reco Jet %i is Et,EtCorr,Eta,Phi (%f ,%f, %f, %f) \n ", RecoJetIdx,recoJet_->et[RecoJetIdx],recoJet_->etCorr[RecoJetIdx] , recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx]);
  for(unsigned int i = 0; i < l1extra_->nCenJets; i++){

    //printf("CenJet Jet is Pt,Eta,Phi ( %f, %f, %f) \n ", l1extra_->cenJetEt[i] , l1extra_->cenJetEta[i], l1extra_->cenJetPhi[i]);

    if( deltaR(recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx], l1extra_->cenJetEta[i], l1extra_->cenJetPhi[i]) < minDeltaR)
    {
      minDeltaR = deltaR(recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx], l1extra_->cenJetEta[i], l1extra_->cenJetPhi[i]);
      //printf( "Delta R = %f \n", minDeltaR);
      if(minDeltaR < .5){
      matchedJet.first = 0;
      matchedJet.second = i;
      }
    }
  }
  //printf("Reco Jet is Pt,Eta,Phi ( %f, %f, %f) \n ", recoJet_->etCorr[RecoJetIdx] , recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx]);
  for(unsigned int i = 0; i < l1extra_->nTauJets; i++){
    //printf("Tau Jet is Pt,Eta,Phi ( %f, %f, %f) \n ", l1extra_->tauJetEt[i] , l1extra_->tauJetEta[i], l1extra_->tauJetPhi[i]);
 if( deltaR(recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx], l1extra_->tauJetEta[i], l1extra_->tauJetPhi[i]) < minDeltaR)
    {
      minDeltaR = deltaR(recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx], l1extra_->tauJetEta[i], l1extra_->tauJetPhi[i]);
      //printf( "Delta R = %f \n", minDeltaR);
      if(minDeltaR < .5){
      matchedJet.first = 1;
      matchedJet.second = i;
      }
    }
  }
  for(unsigned int i = 0; i < l1extra_->nFwdJets; i++){
      //printf("FWD Jet is Pt,Eta,Phi ( %f, %f, %f) \n ", l1extra_->fwdJetEt[i] , l1extra_->fwdJetEta[i], l1extra_->fwdJetPhi[i]);

    if( deltaR(recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx], l1extra_->fwdJetEta[i], l1extra_->fwdJetPhi[i]) <minDeltaR)
    {
      minDeltaR = deltaR(recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx], l1extra_->fwdJetEta[i], l1extra_->fwdJetPhi[i]);
      //printf( "Delta R = %f \n", minDeltaR);
      if(minDeltaR < .5){
      matchedJet.first = 2;
      matchedJet.second = i;
      }
    }
  }
  return matchedJet;
}


std::pair<int,int> L1TowerAnalysis::ReturnMatchedemuJet(int RecoJetIdx){
  double minDeltaR = 9999.;
  std::pair <int,int> matchedJet(-1,-1);
  //printf("Reco Jet %i is Et,EtCorr,Eta,Phi (%f ,%f, %f, %f) \n ", RecoJetIdx,recoJet_->et[RecoJetIdx],recoJet_->etCorr[RecoJetIdx] , recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx]);
  for(unsigned int i = 0; i < l1emuextra_->nCenJets; i++){

    //printf("CenJet Jet is Pt,Eta,Phi ( %f, %f, %f) \n ", l1emuextra_->cenJetEt[i] , l1emuextra_->cenJetEta[i], l1emuextra_->cenJetPhi[i]);

    if( deltaR(recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx], l1emuextra_->cenJetEta[i], l1emuextra_->cenJetPhi[i]) < minDeltaR)
    {
      minDeltaR = deltaR(recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx], l1emuextra_->cenJetEta[i], l1emuextra_->cenJetPhi[i]);
      //printf( "Delta R = %f \n", minDeltaR);
      if(minDeltaR < .5){
      matchedJet.first = 0;
      matchedJet.second = i;
      }
    }
  }
  //printf("Reco Jet is Pt,Eta,Phi ( %f, %f, %f) \n ", recoJet_->etCorr[RecoJetIdx] , recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx]);
  for(unsigned int i = 0; i < l1emuextra_->nTauJets; i++){
    //printf("Tau Jet is Pt,Eta,Phi ( %f, %f, %f) \n ", l1emuextra_->tauJetEt[i] , l1emuextra_->tauJetEta[i], l1emuextra_->tauJetPhi[i]);
 if( deltaR(recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx], l1emuextra_->tauJetEta[i], l1emuextra_->tauJetPhi[i]) < minDeltaR)
    {
      minDeltaR = deltaR(recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx], l1emuextra_->tauJetEta[i], l1emuextra_->tauJetPhi[i]);
      //printf( "Delta R = %f \n", minDeltaR);
      if(minDeltaR < .5){
      matchedJet.first = 1;
      matchedJet.second = i;
      }
    }
  }
  for(unsigned int i = 0; i < l1emuextra_->nFwdJets; i++){
      //printf("FWD Jet is Pt,Eta,Phi ( %f, %f, %f) \n ", l1emuextra_->fwdJetEt[i] , l1emuextra_->fwdJetEta[i], l1emuextra_->fwdJetPhi[i]);

    if( deltaR(recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx], l1emuextra_->fwdJetEta[i], l1emuextra_->fwdJetPhi[i]) <minDeltaR)
    {
      minDeltaR = deltaR(recoJet_->eta[RecoJetIdx], recoJet_->phi[RecoJetIdx], l1emuextra_->fwdJetEta[i], l1emuextra_->fwdJetPhi[i]);
      //printf( "Delta R = %f \n", minDeltaR);
      if(minDeltaR < .5){
      matchedJet.first = 2;
      matchedJet.second = i;
      }
    }
  }
  return matchedJet;
}





double L1TowerAnalysis::ReturnMatchedQuantity( std::pair<int,int> matchJet,int Quantity){
  if(Quantity == 1){
    if(matchJet.first== 0){
      return l1extra_->cenJetEt[matchJet.second];
    }
    if(matchJet.first== 1){
      //cout << " matched jet ET is " << l1extra_->tauJetEt[matchJet.second] << endl;
      return l1extra_->tauJetEt[matchJet.second];
    }
    if(matchJet.first == 2 ){
      return l1extra_->fwdJetEt[matchJet.second];
    }
    if(matchJet.first==-1){
      return false;}
    }

    if(Quantity==2){
        if(matchJet.first==0){
          return l1extra_->cenJetEta[matchJet.second];
        }
        if(matchJet.first==1){
          return l1extra_->tauJetEta[matchJet.second];
        }
        if(matchJet.first==2){
          return l1extra_->fwdJetEta[matchJet.second];
        }
        if(matchJet.first==-1){
          return false;
        }
    }

    if(Quantity==3){
        if(matchJet.first==0){
          return l1extra_->cenJetPhi[matchJet.second];
        }
        if(matchJet.first==1){
          return l1extra_->tauJetPhi[matchJet.second];
        }
        if(matchJet.first==2){
          return l1extra_->fwdJetPhi[matchJet.second];
        }
        if(matchJet.first==-1){
          return false;
        }
    }
    return false;
  }

double L1TowerAnalysis::ReturnMatchedemuQuantity( std::pair<int,int> matchJet,int Quantity){
  if(Quantity == 1){
    if(matchJet.first== 0){
      return l1emuextra_->cenJetEt[matchJet.second];
    }
    if(matchJet.first== 1){
      //cout << " matched jet ET is " << l1emuextra_->tauJetEt[matchJet.second] << endl;
      return l1emuextra_->tauJetEt[matchJet.second];
    }
    if(matchJet.first == 2 ){
      return l1emuextra_->fwdJetEt[matchJet.second];
    }
    if(matchJet.first==-1){
      return false;}
    }

    if(Quantity==2){
        if(matchJet.first==0){
          return l1emuextra_->cenJetEta[matchJet.second];
        }
        if(matchJet.first==1){
          return l1emuextra_->tauJetEta[matchJet.second];
        }
        if(matchJet.first==2){
          return l1emuextra_->fwdJetEta[matchJet.second];
        }
        if(matchJet.first==-1){
          return false;
        }
    }

    if(Quantity==3){
        if(matchJet.first==0){
          return l1emuextra_->cenJetPhi[matchJet.second];
        }
        if(matchJet.first==1){
          return l1emuextra_->tauJetPhi[matchJet.second];
        }
        if(matchJet.first==2){
          return l1emuextra_->fwdJetPhi[matchJet.second];
        }
        if(matchJet.first==-1){
          return false;
        }
    }
    return false;
  }

  bool L1TowerAnalysis::PassHLT( TString TrigBit){
  std::vector<TString>::iterator TrigList = event_->hlt.begin();
  std::vector<TString>::iterator TrigEnd = event_->hlt.end();
  TString ne = TrigBit.Chop();
  for( ; TrigList!=TrigEnd ; ++ TrigList){
    // std::cout << *TrigList << " Compared  " << ne << std::endl;
    if( TrigList->Contains( ne ) ) {
     return true; }
     }
  return false;
 }
  bool L1TowerAnalysis::LooseID( int Jet ){
    if( ( recoJet_->eEMF[Jet]>0.01  ) &&  ( recoJet_->n90hits[Jet] > 1 ) && ( recoJet_->fHPD[Jet] < 0.98 ) ) return true;
    else return false;
  }





       double L1TowerAnalysis::MaxL1Et(){
         double a = 0.;
         double b = 0.;
         double c = 0.;

       if((l1extra_->cenJetEt).size() > 0) a = l1extra_->cenJetEt[0];
       if((l1extra_->tauJetEt).size() > 0) b = l1extra_->tauJetEt[0];
       if((l1extra_->fwdJetEt).size() > 0) c = l1extra_->fwdJetEt[0];
         if(max(a,b) >c){
           return ( a > b ? a : b );
         }
         else return c;
       }

      int L1TowerAnalysis::leadingOfflineJet(){
        double maxEn = 0.;
        int index = -1;

        for(int i = 0; i < (recoJet_->etCorr).size(); ++i)
        {
          if(recoJet_->etCorr[i] > maxEn){
            // //printf("Highest RecoJet is %d , with an Et of %f \n",i,recoJet_->etCorr[i]);
            maxEn = recoJet_->etCorr[i];
            index = i;
          }
        }
        return index;
      }





