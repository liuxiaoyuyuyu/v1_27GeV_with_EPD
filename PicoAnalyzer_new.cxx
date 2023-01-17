#include "PicoAnalyzer_new.h"
#include "StPicoEvent/StPicoEpdHit.h"
#include "StPicoEvent/StPicoBbcHit.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StEpdUtil/StEpdGeom.h"
#include "StEpdUtil/StBbcGeom.h"
#include "StEpdUtil/StEpdEpFinder.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoBTofHit.h"
#include "StRefMultCorr/StRefMultCorr.h"
#include "StRefMultCorr/CentralityMaker.h"


#include "TChain.h"
//#include "TTree.h"
//#include "TLeaf.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TString.h"
#include "TLorentzVector.h"

#include <iostream>
#include <fstream>
using namespace std;

ClassImp(PicoAnalyzer_new)                     //  Macro for CINT compatability


//=================================================
PicoAnalyzer_new::PicoAnalyzer_new(TString FileNameBase):mpTMin(0.15),mpTMax(2.0),mEtaMin(-1.2),mEtaMax(1.2),mNpTbin(12),mNVzbin(16),mNEtabin(24),
mNPhibin(100),mVzMin(-40.0),mVzMax(40.0),mVtxR(1.0),mDiffVzVPD(3.0),mNhitsfit(15),mNhitsfitratio(0.52),mDCAcut(3.0),mFourierOrder(8),
mNTPCSubEvents(2),mNEPDSubEvents(10),mEPDMax(2.0),mEPDthresh(0.3),mpTbound(0.425),mPionSigma(0.012),mKaonSigma(0.012),mProtonSigma(0.012),mEtaMaxv1(0.8),mEtaMinv1(-0.8){
  mFileNameBase = FileNameBase;

  mPicoDst=0;
  mEpdHits=0;
  mBbcHits=0;
  mTracks=0;
  mEventClonesArray=0;

  mRunId=0;
  mRunEt=0;
  mRunCollisionSystem=0;

  mEpdGeom = new StEpdGeom;
  mRan = new TRandom3;
  mRan->GetSeed();
}

//=================================================
PicoAnalyzer_new::~PicoAnalyzer_new(){
  /* no-op */
}


//=================================================
//void PicoAnalyzer_new::SetPicoDst(TTree* PicoDst){
void PicoAnalyzer_new::SetPicoDst(TChain* PicoDst){
  mPicoDst        = PicoDst;

  mEpdHits = new TClonesArray("StPicoEpdHit");
  mBbcHits = new TClonesArray("StPicoBbcHit");
  mTracks  = new TClonesArray("StPicoTrack");
  mEventClonesArray = new TClonesArray("StPicoEvent");
  mTraits = new TClonesArray("StPicoBTofPidTraits");

  mPicoDst->SetBranchStatus("*",0);         // turns OFF all branches  (speeds it up :-)
  unsigned int found;
  mPicoDst->SetBranchStatus("EpdHit*",1,&found);   // note you need the asterisk
  cout << "EpdHit Branch returned found= " << found << endl;
  mPicoDst->SetBranchAddress("EpdHit",&mEpdHits);
  mPicoDst->SetBranchStatus("Event*",1,&found);
  cout << "Event Branch returned found= " << found << endl;
  mPicoDst->SetBranchAddress("Event",&mEventClonesArray);


  mPicoDst->SetBranchStatus("Track*",1,&found);
  cout << "Track Branch returned found= " << found << endl;
  mPicoDst->SetBranchAddress("Track",&mTracks);
  
  mPicoDst->SetBranchStatus("BTofPidTraits*",1,&found);
  cout << "BTofPidTraits Branch returned found= " << found << endl;
  mPicoDst->SetBranchAddress("BTofPidTraits",&mTraits);
}

//=================================================
short PicoAnalyzer_new::Init(char const* TPCWeightFile, char const* TPCShiftFile){

  // ------------------- for the EP finder ------------------
  TString EpFinderOutputName = mFileNameBase;
  EpFinderOutputName += "EpFinderCorrectionsOUTPUT.root";

  mEpFinder = new StEpdEpFinder(144,EpFinderOutputName.Data(),"EPDcorrection.root");// EventType:centrality and Vz, read the INPUT file: EPDcorrection.root 
  mEpFinder->SetnMipThreshold(0.3);
  mEpFinder->SetMaxTileWeight(2.0);
  mEpFinder->SetEpdHitFormat(2);     // 2=pico

  double lin[9] = {-1.950, -1.900, -1.850, -1.706, -1.438, -1.340, -1.045, -0.717, -0.700};
  double cub[9] = {0.1608, 0.1600, 0.1600, 0.1595, 0.1457, 0.1369, 0.1092, 0.0772, 0.0700};
  TH2D wt("Order1etaWeight","Order1etaWeight",100,1.5,6.5,144,0,144);
  for (int ix=1; ix<101; ix++){
    for (int icent=0; icent<9; icent++){
      for(int ivz=0;ivz<16;ivz++){
        double eta = wt.GetXaxis()->GetBinCenter(ix);
        int iy=icent*16+ivz+1;
        wt.SetBinContent(ix,iy,lin[icent]*eta+cub[icent]*pow(eta,3));
      }
    }
  }
  mEpFinder->SetEtaWeights(1,wt);
  cout<<"etaweight set"<<endl;  

  // --------------------------------------------------------

//--------------Prepare the constant for weighting the TPC tracks------------
  mNumberOfTrackTypes =mNpTbin*mNVzbin*mNEtabin*2;
  mNumberOfTrackTypes2 =mNpTbin*mNVzbin*2;
  
  TString OutputRootFileName = mFileNameBase;
  OutputRootFileName += "_FlowHists.root";
  mHistoFile = new TFile(OutputRootFileName,"RECREATE");

  for(int cent=0;cent<9;cent++){
    mVz[cent]= new TH1D(Form("VzCent%d",cent),Form("VzCent%d",cent),16,-40,40);
  }

//---------------Make histograms for checking the flattening of EP---------
//---------------and calculate the resolution------------------------------
    for(int cent=0;cent<9;cent++){
    for(int ivz=0;ivz<16;ivz++){
      for(int iorder=0;iorder<_PsiOrderMax;iorder++){
        mEPDFullPsiWeighted[cent][ivz][iorder] = new TH1F(Form("EPDFullPsi%dWeightedCent%dVz%d",iorder+1,cent,ivz),Form("EPDFullPsi%dWeightedCent%dVz%d",iorder+1,cent,ivz),100,0.0,2.0*TMath::Pi()/(iorder+1.0));
        mEPDFullPsiShifted[cent][ivz][iorder] = new TH1F(Form("EPDFullPsi%dShiftedCent%dVz%d",iorder+1,cent,ivz),Form("EPDFullPsi%dShiftedCent%dVz%d",iorder+1,cent,ivz),100,0.0,2.0*TMath::Pi()/(iorder+1.0));
        mEPDPsiEvsWWeighted[cent][ivz][iorder] = new TH2F(Form("EPDPsi%dEvsWWeightedCent%dVz%d",iorder+1,cent,ivz),Form("EPDPsi%dEvsWWeightedCent%dVz%d",iorder+1,cent,ivz),100,0.0,2.0*TMath::Pi()/(iorder+1.0),100,0.0,2.0*TMath::Pi()/(iorder+1.0));
        mEPDPsiEvsWShifted[cent][ivz][iorder] = new TH2F(Form("EPDPsi%dEvsWShiftedCent%dVz%d",iorder+1,cent,ivz),Form("EPDPsi%dEvsWShiftedCent%dVz%d",iorder+1,cent,ivz),100,0.0,2.0*TMath::Pi()/(iorder+1.0),100,0.0,2.0*TMath::Pi()/(iorder+1.0));
        mResolution[cent][ivz][iorder] = new TProfile(Form("ResolutionPsi%dCent%dVz%d",iorder,cent,ivz),Form("ResolutionPsi%dCent%dVz%d",iorder,cent,ivz),10,0.5,10.5);//x axis corresponds to the <cos()> of different combinations of EP
        mResolution[cent][ivz][iorder]->Sumw2();
        mTPCPsiDisWeighted[cent][ivz][iorder] = new TH1F(Form("TPCPsi%dDisWeightedCent%dVz%d",iorder+1,cent,ivz),Form("TPCPsi%dDisWeightedCent%dVz%d",iorder+1,cent,ivz),100,-TMath::Pi(),TMath::Pi());
        mTPCPsiDisShifted[cent][ivz][iorder] = new TH1F(Form("TPCPsi%dDisShiftedCent%dVz%d",iorder+1,cent,ivz),Form("TPCPsi%dDisShifteddCent%dVz%d",iorder+1,cent,ivz),100,-TMath::Pi(),TMath::Pi());
      }
      for(int ew=0;ew<2;ew++){
        for(int itile=0;itile<24;itile++){
          mThreeDTile[cent][ivz][ew][itile] = new TH3F(Form("dNddeltaphidnMIPdetaCent%dVz%dEW%dTile%d",cent,ivz,ew,itile),Form("dNddeltaphidnMIPdetaCent%dVz%dEW%dTile%d",cent,ivz,ew,itile),24,-1.0*TMath::Pi(),TMath::Pi(),80,0.0,8.0,16,0.5,16.5);//(phi-EP),nMIP,Ring Id 
          mThreeDTile[cent][ivz][ew][itile]->Sumw2();
        }
      }
    }
  }

  for(int cent=0;cent<9;cent++){
    for(int ew=0;ew<2;ew++){
      mAveEta[cent][ew] = new TProfile2D(Form("AveEtaCent%dEW%d",cent,ew),Form("AveEtaCent%dEW%d",cent,ew),16,0.5,16.5,16,0,16);
      mAveEta[cent][ew]->Sumw2();
    }
  }

//---------------Open histograms for getting TPCPhiWeight-----------------
  mTPCWeightFile = new TFile(TPCWeightFile,"READ");

  if(mTPCWeightFile->IsZombie()){
    std::cout<<"TPCWeightFile doesn't exist, didn't apply any phi weight :("<<std::endl;
    for(int icent=0;icent<9;icent++){
      mTPCPhiWeightInput[icent] = 0;
      mTPCPhiAveragedInput[icent] = 0;
      mTPCEtaWeightInput[icent] = 0;
      mTPCEtaAveragedInput[icent] = 0;
    }
  }
  else{
    for(int icent=0;icent<9;icent++){
      mTPCPhiWeightInput[icent] = (TH2D*)mTPCWeightFile->Get(Form("TPCPhiWeightOutputCent%d",icent));
      mTPCPhiAveragedInput[icent] = (TH2D*)mTPCWeightFile->Get(Form("TPCPhiAveragedCent%d",icent));
      mTPCPhiWeightInput[icent]->Divide(mTPCPhiAveragedInput[icent]);
      mTPCEtaWeightInput[icent] = (TH2D*)mTPCWeightFile->Get(Form("TPCEtaWeightOutputCent%d",icent));
      mTPCEtaAveragedInput[icent] = (TH2D*)mTPCWeightFile->Get(Form("TPCEtaAveragedCent%d",icent));
      mTPCEtaWeightInput[icent]->Divide(mTPCEtaAveragedInput[icent]);
    }
  }
//----------------Open histograms for the TPC Psi-shifting----------------
  mTPCShiftFile = new TFile(TPCShiftFile,"READ");

  if(mTPCShiftFile->IsZombie()){
    std::cout<<"TPCShiftFile doesn't exist, didn't do Psi-shifting :("<<std::endl;
    for(int ivz=0;ivz<16;ivz++){
      mTPCCosShift_Input[ivz] = 0;
      mTPCSinShift_Input[ivz] = 0;
    }
  }
  else{
    for(int ivz=0;ivz<16;ivz++){
      mTPCCosShift_Input[ivz] = (TProfile3D*)mTPCShiftFile->Get(Form("TPCCosShiftVz%d",ivz));
      mTPCSinShift_Input[ivz] = (TProfile3D*)mTPCShiftFile->Get(Form("TPCSinShiftVz%d",ivz));
    }
  }

  cout << "Init done\n";
  return 0;
}

//==============================================
// some things might need resetting when there is a new run
void PicoAnalyzer_new::NewRun(int runId){
  //mRunCollisionSystem = WhichSystem(runId); //tend to cause error.
  mRunId = runId;
  //mRunEt += 1; 
}

//=================================================
short PicoAnalyzer_new::Make(int iEvent){
  //----------------- get data --------------
  mPicoDst->GetEntry(iEvent);
  StPicoEvent* event = (StPicoEvent*)((*mEventClonesArray)[0]);
  if (event->runId()!=mRunId){
    NewRun(event->runId());        // some things should be reset when there is a new run loaded
    cout << "New run detected: " << mRunId << " and it is collision system #" << mRunCollisionSystem << endl;
    //cout<<"RunEntry"<<mRunEt<<endl;  
  }

  //----- done getting data; have fun! ------
  if(!(event->isTrigger(610001)||event->isTrigger(610011)||event->isTrigger(610021)||event->isTrigger(610031)||event->isTrigger(610041)||event->isTrigger(610051))) return 0;
  
  //  StThreeVectorF primaryVertex = event->primaryVertex();
  //  TVector3 PV(primaryVertex.x(),primaryVertex.y(),primaryVertex.z());
  TVector3 PV = event->primaryVertex();
  int mRunId = event->runId();
  int RUNYear = 2018;
  float RUNEnergy = 27.;
  int RunDay = floor( (mRunId - (RUNYear-2000)*pow(10,6))/pow(10,3) );
  int DayBinId = RunDay-89;

  double pi = TMath::Pi();

  //-----------------Get the multiplicity by the StRefMulCorr class--------------
  int CentId=-1;
  int mRefMult=event->refMult();
  
  bool ISRefMultCorrBadRun=false;
  mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
  //mRefMultCorr = new StRefMultCorr();
  mRefMultCorr->init(mRunId);
  if(mRefMultCorr->getBeginRun(RUNEnergy,RUNYear)==-1) return 0;
  ISRefMultCorrBadRun=mRefMultCorr->isBadRun(mRunId);
  if(ISRefMultCorrBadRun) return 0;
  mRefMultCorr->initEvent(mRefMult,PV.Z(),event->ZDCx());
  mRefMult = mRefMultCorr->getRefMultCorr();
  CentId = mRefMultCorr->getCentralityBin9();//An integer between 0 (70-80%) and 8 (0-5%)
  //int CentIdmy = FindCent(event->refMult());   // returns an integer between 0 (70-80%) and 8 (0-5%)
  //-------------remove the pile-up events-----------------
  if(mRefMultCorr->passnTofMatchRefmultCut(1.*event->refMult(), 1.*event->nBTOFMatch())!=1) return 0;
  if (fabs(PV.Z())>=mVzMax) return 0;
  if (sqrt(pow(PV.X(),2)+pow(PV.Y(),2))>mVtxR) return 0;
  //if (abs(PV.Z()-event->vzVpd())>mDiffVzVPD) return 0;//Get rid of the VPD cut, Prithwish said it does no good at low energy 06/18/2020
  if (CentId<0) return 0;            // 80-100% - very peripheral


  int VzBin;
  double VzArr[17]={-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40};
  for(int i=0;i<16;i++){
    if(PV.Z()>VzArr[i]&&PV.Z()<=VzArr[i+1]){
      VzBin=i;
      break;
    }
  }
  //Get BinId for StEpdEpFinder
  int EpdId=CentId*16+VzBin;

  StEpdEpInfo result = mEpFinder->Results(mEpdHits,PV,EpdId);  // and now you have all the EP info you could ever want :-)
  //StEpdEpInfo result = mEpFinder->Results(mEpdHits,PV,1);  // testing. and now you have all the EP info you could ever want :-)
  double EpAngle[_PsiOrderMax][3];//EP order,east/west/full
  for(int iorder=0;iorder<_PsiOrderMax;iorder++){
    EpAngle[iorder][0]=result.EastPhiWeightedAndShiftedPsi(iorder+1);
    EpAngle[iorder][1]=result.WestPhiWeightedAndShiftedPsi(iorder+1);
    EpAngle[iorder][2]=result.FullPhiWeightedAndShiftedPsi(iorder+1);
  }

  //mHisto2D[0]->Fill(EpAngle[1],EpAngle[0]);//East vs. West

  mVz[CentId]->Fill(PV.Z());
  //-----------Check the flattening of Psi_EPDFull----------
  for(int iorder=0;iorder<_PsiOrderMax;iorder++){
    mEPDFullPsiWeighted[CentId][VzBin][iorder]->Fill(result.FullPhiWeightedPsi(iorder+1));
    mEPDFullPsiShifted[CentId][VzBin][iorder]->Fill(result.FullPhiWeightedAndShiftedPsi(iorder+1));
    mEPDPsiEvsWWeighted[CentId][VzBin][iorder]->Fill(result.EastPhiWeightedPsi(iorder+1),result.WestPhiWeightedPsi(iorder+1));
    mEPDPsiEvsWShifted[CentId][VzBin][iorder]->Fill(result.EastPhiWeightedAndShiftedPsi(iorder+1),result.WestPhiWeightedAndShiftedPsi(iorder+1));
  }
  
  //------------Prepare Qvectors for TPC EP (used in v1EPD)--------------
  TH1D* PCosPhi;
  TH1D* PSinPhi;
  //PCosPhi = new TH1D(Form("PCosPhi"),Form("PCosPhi"),3,0.5,3.5);
  //PSinPhi = new TH1D(Form("PSinPhi"),Form("PSinPhi"),3,0.5,3.5);//x is three types of Psi_TPC
  PCosPhi = new TH1D(Form("PCosPhi"),Form("PCosPhi"),_PsiOrderMax,0.5,_PsiOrderMax+0.5);
  PSinPhi = new TH1D(Form("PSinPhi"),Form("PSinPhi"),_PsiOrderMax,0.5,_PsiOrderMax+0.5);

  //------------Begin loop over TPC tracks--------------------------
  for(int itrk=0; itrk<mTracks->GetEntries(); itrk++){
    StPicoTrack* track = (StPicoTrack*)((*mTracks)[itrk]);
    if ((track->nHitsFit()<mNhitsfit)||(!track->isPrimary())) continue; // I should check DCA too, but am confused how
    double nHitsFitRatio = track->nHitsFit()*1.0/track->nHitsMax();
    //cout<<nHitsFitRatio<<endl;
    //if (nHitsFitRatio<mNhitsfitratio) continue;//get rid of the nhitsfitratio cut, Prithwish said it is a very old cut used by STAR. 06/18/2020
    TVector3 pMom = track->pMom();//track->gMom() if I want to look at the global tracks.
    double Tphi = pMom.Phi();
    double Teta = pMom.Eta();
    double TPt = pMom.Pt();//
    double dca = track->gDCA(PV).Mag();
    if (dca>mDCAcut) continue;
    //mHisto1D[3]->Fill(TPt);
    int Tch=track->charge();
    double rig=Tch*pMom.Mag();
    double dEdx=track->dEdx();

    if(TMath::Abs(Teta)>=mEtaMaxv1) continue;
    if(TPt<=mpTMin||TPt>=mpTMax) continue;//pT range [0.15,2.0] for Psi_TPC
    //----------Prepare Q to calculate the Psi_TPC-----------------
    int TtrId=FindTrackId(track->charge(),PV.Z(),Teta,TPt);
    int TtrId2=FindTrackId2(track->charge(),PV.Z(),TPt);
    //-------TPC track weighting histos---------
    if(TtrId<=0||TtrId>mNumberOfTrackTypes){
      cout<<"Ah-oh, invalid TrackId :("<<endl;
    }
    if(TtrId2<=0||TtrId2>mNumberOfTrackTypes2){
      cout<<"Ah-oh, invalid TrackId2 :("<<endl;
    }
    
    double TrackPhiWeight=1.0;
    if(mTPCPhiWeightInput[CentId]!=0) TrackPhiWeight=1.0/mTPCPhiWeightInput[CentId]->GetBinContent(mTPCPhiWeightInput[CentId]->FindBin(Tphi,TtrId));
    if(!TMath::Finite(TrackPhiWeight)) continue;
    if(TrackPhiWeight>10.0) continue;
    double TrackEtaWeight=1.0;
    if(mTPCEtaWeightInput[CentId]!=0) TrackEtaWeight=1.0/mTPCEtaWeightInput[CentId]->GetBinContent(mTPCEtaWeightInput[CentId]->FindBin(Teta,TtrId2));
    if(!TMath::Finite(TrackEtaWeight)) continue; 
    double TrackWeight=TrackPhiWeight*TrackEtaWeight;  
    if(!TMath::Finite(TrackWeight)) continue;
    
    double TPCWeight[_PsiOrderMax];
    TPCWeight[0] = TrackWeight*(-1.0*Teta);//weight the TPC tracks with (-eta) for Psi1
    //TPCWeight[1] = TrackPhiWeight*TPt;//weight the TPC tracks with pT for Psi2

    for(int iorder=0;iorder<_PsiOrderMax;iorder++){
      PCosPhi->Fill(iorder+1,TPCWeight[iorder]*cos((double)(iorder+1)*Tphi));
      PSinPhi->Fill(iorder+1,TPCWeight[iorder]*sin((double)(iorder+1)*Tphi));
    }

  }  
  //-------------End loop over TPC tracks-----------------------

  //-------------Now let's play with the EP---------------------
  //double TPCPsi1[3]={0.0};// Three types of Psi_TPC: 0_Full, 1_pos, 2_neg
  double TPCPsi[_PsiOrderMax];
  for(int i=1; i<=_PsiOrderMax;i++){
    TPCPsi[i-1]=-999;
    if(fabs(PSinPhi->GetBinContent(i))>1e-6&&fabs(PCosPhi->GetBinContent(i))>1e-6){
      TPCPsi[i-1] = atan2(PSinPhi->GetBinContent(i),PCosPhi->GetBinContent(i))/(double)i;
    }
    //atan2 returns angle between -pi and pi, remember to divide by i!!!! 06/16/2020
  }
  
  //-------------Shift the TPC EP if the shifting histos exist------------
  double TPCPsiShifted[_PsiOrderMax];
  for(int i=0;i<_PsiOrderMax;i++){
    TPCPsiShifted[i]=TPCPsi[i];
  }

  int TPCPsiId=11-int(mEtaMaxv1/0.1);

  if(mTPCCosShift_Input[VzBin]!=0&&mTPCSinShift_Input[VzBin]!=0){
    double deltapsi[_PsiOrderMax]={0.0};
    for(int i=1;i<=_PsiOrderMax;i++){
      for(int k=1;k<=mFourierOrder;k++){
        double tmp=(double)(i*k);
        deltapsi[i-1]+=2.0*(-mTPCSinShift_Input[VzBin]->GetBinContent(CentId+1,k,TPCPsiId)*cos(tmp*TPCPsi[i-1])+mTPCCosShift_Input[VzBin]->GetBinContent(CentId+1,k,TPCPsiId)*sin(tmp*TPCPsi[i-1]))/tmp;
        //z=1,[-1.0,1.0];z=2,[-0.9,0.9];z=3,[-0.8,0.8]....
      }
      if(TPCPsiShifted[i-1]>-100){
        TPCPsiShifted[i-1]+=deltapsi[i-1];
        if(TPCPsiShifted[i-1]<-pi/(double)i) TPCPsiShifted[i-1]+=2*pi/(double)i;
        if(TPCPsiShifted[i-1]>pi/(double)i) TPCPsiShifted[i-1]-=2*pi/(double)i;
      }
    }  
  }

  //-------------Fill the TProfiles for calculating the Resolution---------
  for(int i=0;i<_PsiOrderMax;i++){
    mResolution[CentId][VzBin][i]->Fill(1,cos((double)(i+1)*(EpAngle[i][0]-TPCPsiShifted[i])));
    mResolution[CentId][VzBin][i]->Fill(2,cos((double)(i+1)*(EpAngle[i][1]-TPCPsiShifted[i])));
    mResolution[CentId][VzBin][i]->Fill(3,cos((double)(i+1)*(EpAngle[i][0]-EpAngle[i][1])));
  }

  //----Fill the histograms for checking the flattening of the TPC EPs-----
  for(int i=0;i<_PsiOrderMax;i++){
    mTPCPsiDisWeighted[CentId][VzBin][i]->Fill(TPCPsi[i]);
    mTPCPsiDisShifted[CentId][VzBin][i]->Fill(TPCPsiShifted[i]);
  }

  //--------------Begin loop over EPD hits---------------------------------
  for(int hit=0;hit<mEpdHits->GetEntries();hit++){
    int tileId, ring, TT, PP, EW, ADC;
    float nMip;
    StPicoEpdHit* epdHit = (StPicoEpdHit*)((*mEpdHits)[hit]);
    tileId=epdHit->id();
    EW = (tileId<0)?0:1;
    ring = epdHit->row();//[1,16]
    TT = epdHit->tile();
    PP = epdHit->position();
    ADC = epdHit->adc();
    nMip = epdHit->nMIP();

    int TTxyId;
    if(EW==0){
      int OE=(TT%2)?2:1;//odd 2, even 1
      int PPId=(PP>9)?(21-PP):(9-PP);
      TTxyId=PPId*2+OE;
    }
    else{
      int OE=(TT%2)?1:2;//odd 1, even 2
      int PPId=(PP<4)?(PP+8):(PP-4);
      TTxyId=PPId*2+OE;
    }
    
    //if (nMip<mEPDthresh) continue;
    //double TileWeight = (nMip<3.0)?nMip:3.0;//Note: here a different nMIPMax from the EpFinder was used
    TVector3 StraightLine = mEpdGeom->RandomPointOnTile(tileId) - PV;
    double Hphi = StraightLine.Phi();
    double Heta = StraightLine.Eta();

    mAveEta[CentId][EW]->Fill(ring,VzBin,Heta);

    double deltaphi[_PsiOrderMax];//order, reference: 0-TPC,1-the other side of EPD
    for(int iorder=0;iorder<_PsiOrderMax;iorder++){
      deltaphi[iorder]=Hphi-TPCPsiShifted[iorder];
    }

    for(int iorder=0;iorder<_PsiOrderMax;iorder++){
      double tmpn=pi/(double)(iorder+1);
      while(deltaphi[iorder]<-tmpn||deltaphi[iorder]>tmpn){
        if(deltaphi[iorder]<-tmpn) deltaphi[iorder]+=2.0*tmpn;
        else if(deltaphi[iorder]>tmpn) deltaphi[iorder]-=2.0*tmpn;
      }
    }
    
    mThreeDTile[CentId][VzBin][EW][TTxyId-1]->Fill(deltaphi[0],nMip,ring); 

  }
  //-----------------End looping over EPD hits------------------------

  delete PCosPhi;
  delete PSinPhi;

  return 0;
}
//=================================================
short PicoAnalyzer_new::Finish(){

  cout << "In PicoAnalyzer_new::Finish - calling StEpdEpFinder::Finish()\n";
  mEpFinder->Finish();

  cout << "I have called it\n";

  mTPCWeightFile->Close();
  mTPCShiftFile->Close();
  cout<<"Hello, I just want some weighting hisots"<<endl;

  mHistoFile->Write();
  mHistoFile->Close();

  cout << "Finish!!\n\n";

  return 0;
}
/*
bool PicoAnalyzer_new::Runlist(int runId){
    bool runpass=0;
    int runlist[]={19130071,19130077,19130078,19130079,19130084,19130085,19130086,19131001,19131003,19131004,19131005,19131006,19131007,19131009,19131010,19131012,19131013,19131014,19131015,19131016,19131019,19131020,19131026,19131027,19131030,19131031,19131037,19131039,19131040,19131041,19131042,19131044,19131045,19131046,19131048,19131049,19131050,19131051,19131052,19131054,19131055,19131056,19131057,19131060,19131061,19131062,19132001,19132004,19132005,19132006,19132007,19132008,19132011,19132013,19132014,19132016,19132017,19132019,19132020,19132021,19132022,19132026,19132027,19132029,19132030,19132031,19132032,19132034,19132035,19132036,19132037,19132038,19132039,19132043,19132044,19132045,19132046,19132047,19132048,19132063,19132064,19132065,19132070,19132071,19132073,19132074,19132075,19132076,19132078,19132079,19132080,19132081,19132082,19132083,19133002,19133003,19133004,19133005,19133008,19133009,19133010,19133012,19133013,19133014,19133016,19133017,19133018,19133021,19133022,19133023,19133025,19133026,19133027,19133028,19133030,19133031,19133032,19133033,19133038,19133039,19133040,19133041,19133043,19133044,19133045,19133047,19133048,19133049,19133050,19133053,19133054,19133055,19133056,19133058,19133059,19133060,19133061,19134002,19134003,19134004,19134005,19134007,19134008,19134009,19134010,19134011,19134013,19134014,19134016,19134017,19134018,19134019,19134023,19134024,19134025,19134027,19134028,19134029,19134030,19134035,19134036,19134037,19134038,19134040,19134041,19134042,19134044,19134045,19134046,19134047,19134049,19134050,19135001,19135003,19135004,19135005,19135012,19135013,19135014,19135016,19135017,19135018,19135020,19135021,19135022,19135027,19135028,19135029,19135037,19135038,19135039,19135040,19135042,19135043,19136001,19136003,19136004,19136005,19136007,19136008,19136009,19136011,19136012,19136013,19136014,19136016,19136017,19136018,19136040,19136041,19136042,19136044,19136045,19136046,19136047,19136049,19137001,19137002,19137003,19137004,19137006,19137007,19137008,19137009,19137010,19137011,19137013,19137014,19137015,19137016,19137019,19137020,19137022,19137023,19137024,19137025,19137026,19137027,19137028,19137029,19137036,19137037,19137038,19137040,19137041,19137045,19137047,19137048,19137050,19137051,19137052,19137053,19137054,19137056,19137057,19137058,19137059,19138002,19138003,19138004,19138005,19138006,19138008,19138009,19138010,19138011,19138012,19138014,19138015,19138016,19138019,19138020,19138021,19138025,19138026,19138027,19138028,19139023,19139024,19139026,19139027,19139028,19139032,19139033,19139034,19139037,19139038,19139039,19139041,19139042,19139043,19139044,19139050,19139058,19139063,19139064,19139066,19139067,19139068,19139069,19139071,19139072,19139073,19140002,19140003,19140004,19140006,19140007,19140008,19140010,19140011,19140012,19140014,19140015,19140016,19140017,19140020,19140021,19140022,19140025,19140030,19140031,19140035,19140036,19140037,19140040,19140041,19140042,19140043,19140045,19140046,19140047,19140051,19140052,19140053,19140055,19140056,19141001,19141003,19141004,19141005,19141006,19141008,19141009,19141010,19141011,19141013,19141014,19141015,19141018,19141019,19141020,19141021,19141022,19141024,19141025,19141026,19141028,19141029,19141030,19141047,19141048,19141049,19141051,19141052,19141053,19142001,19142002,19142003,19142005,19142006,19142007,19142008,19142010,19142011,19142012,19142014,19142015,19142016,19142018,19142019,19142020,19142021,19142022,19142027,19142034,19142035,19142038,19142039,19142041,19142045,19142048,19142049,19142053,19142054,19142055,19142057,19142058,19142059,19142062,19142064,19142065,19142068,19143001,19143003,19143006,19143007,19143008,19143009,19143010,19143011,19143012,19143013,19143014,19143015,19143016,19143017,19144012,19144013,19144014,19144018,19144019,19144020,19144024,19144025,19144026,19144031,19144032,19144033,19144036,19144037,19144038,19144042,19144043,19144044,19144046,19144047,19145001,19145004,19145005,19145006,19145008,19145009,19145010,19145011,19145013,19145014,19145015,19145017,19145019,19145020,19145028,19145031,19145034,19145035,19145036,19145038,19145039,19145040,19145042,19145043,19145044,19145047,19145048,19145050,19146002,19146003,19146004,19146006,19146007,19146008,19146009,19146012,19146013,19146014,19146016,19146017,19146019,19146020,19146024,19146025,19146026,19147007,19147008,19147009,19147014,19147015,19147016,19147021,19147022,19147023,19147025,19147026,19147027,19147029,19147030,19147031,19147033,19147034,19147035,19147038,19147039,19147040,19147042,19147043,19147044,19147046,19147047,19147048,19148002,19148003,19148004,19148007,19148008,19148009,19148011,19148012,19148013,19148015,19148016,19148017,19148021,19148022,19148023,19148024,19148050,19148051,19148052,19149002,19149003,19149007,19149008,19149009,19149012,19149013,19149014,19149015,19149017,19149018,19149023,19149024,19149025,19149030,19149031,19149032,19149035,19149037,19149038,19149039,19149042,19149043,19149044,19149046,19149047,19149048,19149050,19149051,19149052,19149054,19149055,19150001,19150005,19150006,19150007,19150009,19150010,19150011,19150013,19150014,19150016,19150017,19150018,19155057,19155058,19156001,19156002,19156004,19156005,19156006,19156008,19156009,19156011,19156013,19156014,19156015,19156018,19156019,19156030,19156031,19156032,19156033,19156042,19156043,19156044,19156045,19156046,19156047,19157002,19157003,19157004,19157006,19157007,19157008,19157012,19157013,19157015,19157017,19157018,19158020,19158059,19158060,19158062,19159001,19159003,19159004,19159006,19159007,19159008,19159010,19159011,19159012,19159014,19159015,19159016,19159018,19159019,19159021,19159023,19159024,19159025,19159039,19159040,19159041,19160002,19160003,19160004,19160006,19160007,19160008,19160010,19160011,19160012,19160014,19160015,19160016,19160019,19160020,19160021,19160023,19160024,19160025,19160027,19160028,19160029,19161003,19161004,19161005,19161007,19161008,19161009,19161011,19161012,19161013,19161015,19161016,19161017,19161048,19161049,19161050,19161051,19161053,19162002,19162006,19162008,19162009,19162010,19162013,19162014,19162015,19162017,19162018,19162019,19162021,19162022,19162023,19162030,19162031,19162032,19162036,19162037,19162038,19162040,19162041,19162042,19163006,19163007,19163008,19163010,19163011,19163012,19163014,19163015,19163016,19163018,19163019,19163020,19163039,19163040,19163041,19164004,19164005,19164006,19164008,19164009,19164010,19164011,19164012,19164014,19164015,19164016,19164017,19165002,19165003,19165004,19165006,19165007,19165008,19165009,19165011,19165012,19165013,19165015,19165016,19165017,19165018,19165020,19165021,19165027,19165028,19165029,19166006,19166007,19166008,19166010,19166011,19166012,19166014,19166015,19167003,19167004,19167005,19167007,19167008,19167009,19167011,19167012,19167013,19167015,19167016,19167017,19167021,19167022,19167023,19167025,19167026,19167027,19167028,19167030,19167031,19167032,19167034,19167035,19167042,19167044,19167045,19167046,19167049,19168025,19168026,19168028,19168029,19168030,19168033,19168034,19168036,19168038,19168039,19168040};
    for(int i=0;i<788;i++){
      if(runId==runlist[i]){
        runpass=1;
        break;
      }
    }
    return runpass;
}
*/

//TrackId starts at 1!
int PicoAnalyzer_new::FindTrackId(int Trkch,double TrkVz,double TrkEta,double TrkpT){
  int TrackId=0;
  int chId=(Trkch==-1)?0:1;
  int VzId=-1;
  if(TrkVz<0.0) VzId=-(int)(mNVzbin*fabs(TrkVz)/(mVzMax-mVzMin))+mNVzbin/2-1;
  else VzId=(int)(mNVzbin*fabs(TrkVz)/(mVzMax-mVzMin))+mNVzbin/2;
  int EtaId=-1;
  if(TrkEta<0.0) EtaId=-(int)(mNEtabin*fabs(TrkEta)/(mEtaMax-mEtaMin))+mNEtabin/2-1;
  else EtaId=(int)(mNEtabin*fabs(TrkEta)/(mEtaMax-mEtaMin))+mNEtabin/2;
  int pTId=-1;
  if(TrkpT>1.25) pTId=11;
  else pTId=(int)(11*fabs(TrkpT-mpTMin)/(1.25-mpTMin));

  TrackId=chId*mNpTbin*mNVzbin*mNEtabin+VzId*mNpTbin*mNEtabin+EtaId*mNpTbin+pTId+1;
  return TrackId;
}

int PicoAnalyzer_new::FindTrackId2(int Trkch,double TrkVz,double TrkpT){
  int TrackId=0;
  int chId=(Trkch==-1)?0:1;
  int VzId=-1;
  if(TrkVz<0.0) VzId=-(int)(mNVzbin*fabs(TrkVz)/(mVzMax-mVzMin))+mNVzbin/2-1;
  else VzId=(int)(mNVzbin*fabs(TrkVz)/(mVzMax-mVzMin))+mNVzbin/2;
  int pTId=-1;
  if(TrkpT>1.25) pTId=11;
  else pTId=(int)(11*fabs(TrkpT-mpTMin)/(1.25-mpTMin));

  TrackId=chId*mNpTbin*mNVzbin+VzId*mNpTbin+pTId+1;

  return TrackId;
}




