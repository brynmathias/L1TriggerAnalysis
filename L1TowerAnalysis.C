#include "./L1TowerAnalysis.hh"
// --------------------------------------------------------------------
//                       L1TowerAnalysis macro definition
//
//    macro automatically generated by GenerateMacro.sh script
//    author           : bm409
//    creation date    : Sat Oct 23 14:00:25 CEST 2010
//    last update date :
//    description      :
//
// --------------------------------------------------------------------


// --------------------------------------------------------------------
//                             run function
// --------------------------------------------------------------------

void L1TowerAnalysis::BookHistos() {

  RefJets = new TH1F("RefJet", "RefJetEt",1000,0.,1000);
  ET10   = new TH1F( "JetEt10" , "JetEt",1000, 0., 1000.);
  l1emuTauEt = new TH1F("L1TauJetEt", "tauJet",260,0.,260.);
  L1HT = new TH1F("L1HT", "HT",1000,0.,1000.);
  L1MET = new TH1F("L1MET","MET",1000,0.,1000.);
  AllL1Jets = new TH1F("AllL1Jets", "AllL1Jets",260,0.,260.);
  centralJets = new TH1F("centralJets", "centralJets",260,0.,260.);
  NonForward =  new TH1F("NonForward", "NonForward",260,0.,260.);
  forward = new TH1F("forward", "forward",260,0.,260.);
  QuadJet20  = new TH1F("QuadJet20", "QuadJet20",1,0.,1.);
  l1TauEt = new TH1F("l1TauEt", "QuadJet20",1,0.,1.);
  RefJetsemu = new TH1F("RefJetemu", "RefJetEt",1000,0.,1000);
  ET10emu   = new TH1F( "JetEt10emu" , "JetEt",1000, 0., 1000.);
  l1emuTauEtemu = new TH1F("L1TauJetEtemu", "tauJet",260,0.,260.);
  L1HTemu = new TH1F("L1HTemu", "HT",1000,0.,1000.);
  L1METemu = new TH1F("L1METemu","MET",1000,0.,1000.);
  AllL1Jetsemu = new TH1F("AllL1Jetsemu", "AllL1Jets",260,0.,260.);
  centralJetsemu = new TH1F("centralJetsemu", "centralJets",260,0.,260.);
  forwardemu = new TH1F("forwardemu", "forward",260,0.,260.);
  QuadJet20emu  = new TH1F("QuadJet20emu", "QuadJet20",2,0.,2.);
  l1TauEtemu = new TH1F("l1TauEtemu", "l1TauEtemu",1,0.,1.);
  events = new TH1F("events", "events",2,0.,2.);
  emu_extraJetCorrelation = new TH2F("emu_extraJetCorrelation","emu_extraJetCorrelation;l1Extra;L1EmuExtra",26,0.,260.,26,0.,260.);

  }

// --------------------------------------------------------------------
//                             run function
// --------------------------------------------------------------------
  void L1TowerAnalysis::run(Long64_t nevents, TString outputname, TString TriggerBit) //To run use m.run(events,"SomeoutputName","HLT TriggerName") (leave HLT name as white space to select no trigger)
  {
    // Make output File
    TFile *theFile = new TFile(outputname, "RECREATE");
    theFile->cd();
    //Book Histos

    BookHistos();
  
    //number of events to process
    if (nevents==-1 || nevents>GetEntries()) nevents=GetEntries();
    std::cout << nevents  << " to process ..." << std::endl;
    //loop over the events
    for (Long64_t i= 0; i < nevents; i++)
    {
      //load the i-th event
      Long64_t ientry = LoadTree(i); if (ientry < 0) break;
      GetEntry(i);
      //process progress
      if( i!=0 && ( i%10000 ) ==0 ) {std::cout << "- processing event " << i << "\r" << std::flush;}
      double wgt =1.0;
      if (!TriggerBit.IsWhitespace() && !PassHLT(TriggerBit) ){ continue;}
      events->Fill(1,wgt);

      
      double HT = 0.;
      int NJetsGreater20GeVemu = 0;
      int NJetsGreater20GeV = 0;
      // Check that we have jets in the event.
      if((recoJet_->etCorr).size() < 1) continue;
      L1METemu->Fill(l1emuextra_->met,wgt);

    for( unsigned jet = 0; jet < recoJet_->etCorr.size(); ++jet){
      if( recoJet_->etCorr[jet] < 30.) continue;
      if(!LooseID(jet)) continue;
      if(!MatchemuJet(jet)) continue;
      std::pair<int,int> matchedJet = ReturnMatchedemuJet(jet);

      if( recoJet_->etCorr[jet] > 40.&& fabs(recoJet_->eta[jet]) < 3.) HT+=recoJet_->etCorr[jet];


        
        if( matchedJet.first == 0){
      // for( unsigned k = 0; k < l1emuextra_->cenJetEt.size(); ++k){
          AllL1Jetsemu->Fill(l1emuextra_->cenJetEt[matchedJet.second],wgt);
          centralJetsemu->Fill(l1emuextra_->cenJetEt[matchedJet.second],wgt);
          if( l1emuextra_->cenJetEt[matchedJet.second] > 38. ){NJetsGreater20GeVemu++;}
          NonForward->Fill(l1emuextra_->cenJetEt[matchedJet.second],wgt);          
        } 
        
        
        if(matchedJet.first == 1){
      // for( unsigned k = 0; k < l1emuextra_->tauJetEt.size(); ++k){
          l1emuTauEtemu->Fill(l1emuextra_->tauJetEt[matchedJet.second],wgt);
          AllL1Jetsemu->Fill(l1emuextra_->tauJetEt[matchedJet.second],wgt);
          NonForward->Fill(l1emuextra_->tauJetEt[matchedJet.second],wgt);
          if( l1emuextra_->tauJetEt[matchedJet.second] > 38. ){NJetsGreater20GeVemu++;}
        }

        if(matchedJet.first == 2){
      // for( unsigned k = 0; k < l1emuextra_->fwdJetEt.size(); ++k){
          AllL1Jetsemu->Fill(l1emuextra_->fwdJetEt[matchedJet.second],wgt);
          forwardemu->Fill(l1emuextra_->fwdJetEt[matchedJet.second],wgt);
        }
      }
        if( NJetsGreater20GeVemu > 3 ){
          QuadJet20emu->Fill(1.,wgt);
        }
      // std::cout << "NonEmu stuff" << std::endl;
        L1HT->Fill(l1extra_->ht,wgt); 
        L1MET->Fill(l1extra_->met,wgt);
      //   for( unsigned jet = 0; jet < recoJet_->etCorr.size(); ++jet){
      //     if( recoJet_->etCorr[jet] < 30.) continue;
      //     if(!LooseID(jet)) continue;
      //     if(!MatchJet(jet)) continue;
      //       std::pair<int,int> matchedJet = ReturnMatchedJet(jet);
      //   
      //    
      //       if(matchedJet.first == 1){
      //     l1TauEt->Fill(l1extra_->tauJetEt[k],wgt);
      //     AllL1Jets->Fill(l1extra_->tauJetEt[k],wgt);
      //     if( l1extra_->tauJetEt[k] > 38. ){NJetsGreater20GeV++;}
      //   }
      //   if( matchedJet.first == 0){
      //     AllL1Jets->Fill(l1extra_->cenJetEt[k],wgt);
      //     centralJets->Fill(l1extra_->cenJetEt[k],wgt);
      //     if( l1extra_->cenJetEt[k] > 38. ){NJetsGreater20GeV++;}
      //   }
      //   if(matchedJet.first == 2){
      //     AllL1Jets->Fill(l1extra_->fwdJetEt[k],wgt);
      //     forward->Fill(l1extra_->fwdJetEt[k],wgt);
      //   }
      // }
      //   if( NJetsGreater20GeV > 3 ){
      //     QuadJet20->Fill(1.,wgt);
      //   }

        if(HT > 100.)  L1HTemu->Fill(l1emuextra_->ht,wgt); 
        unsigned nJets = l1extra_->cenJetEt.size();
        if( nJets > l1emuextra_->cenJetEt.size()) nJets = l1emuextra_->cenJetEt.size();
        for( unsigned jet = 0; jet < nJets; ++jet){
        emu_extraJetCorrelation->Fill(l1extra_->cenJetEt[jet],l1emuextra_->cenJetEt[jet],1.);
        }
      }
    // Write and close the file!
      theFile->Write();
      theFile->Close();


    }
