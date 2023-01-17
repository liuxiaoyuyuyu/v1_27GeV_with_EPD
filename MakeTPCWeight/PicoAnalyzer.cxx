//#define _OutputTPCWeightFile
//#define _OutputTPCShiftFile
#define _OutputTPCQAFile

#include "FlowAnalysis/PicoAnalyzer.h"
#include "StPicoEvent/StPicoEpdHit.h"
#include "StPicoEvent/StPicoBbcHit.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StEpdUtil/StEpdGeom.h"
#include "StEpdUtil/StBbcGeom.h"
#include "StEpdUtil/StEpdEpFinder.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"


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

#include <iostream>
#include <fstream>
using namespace std;

ClassImp(PicoAnalyzer)                     //  Macro for CINT compatability


//=================================================
PicoAnalyzer::PicoAnalyzer(TString FileNameBase):mpTMin(0.15),mpTMax(2.0),mEtaMin(-1.2),mEtaMax(1.2),mNpTbin(12),mNVzbin(16),mNEtabin(24),
mNPhibin(100),mVzMin(-40.0),mVzMax(40.0),mVtxR(1.0),mDiffVzVPD(3.0),mNhitsfit(15),mNhitsfitratio(0.52),mDCAcut(3.0),mFourierOrder(8),mEtaMaxAly(1.0){
//"..Eta.." are the thresholds of the weighting histos. "..EtaAly.." are the thresholds for the analysis, e.g EP from TPC 

  mFileNameBase = FileNameBase;

  mPicoDst=0;
  mTracks=0;
  mEventClonesArray=0;

  mRunId=0;
  mRunEt=0;
  mRunCollisionSystem=0;

  mRan = new TRandom3;
  mRan->GetSeed();
}

//=================================================
PicoAnalyzer::~PicoAnalyzer(){
  /* no-op */
}


//=================================================
//void PicoAnalyzer::SetPicoDst(TTree* PicoDst){
void PicoAnalyzer::SetPicoDst(TChain* PicoDst){
  mPicoDst        = PicoDst;

  //mEpdHits = new TClonesArray("StPicoEpdHit");
  //mBbcHits = new TClonesArray("StPicoBbcHit");
  mTracks  = new TClonesArray("StPicoTrack");
  mEventClonesArray = new TClonesArray("StPicoEvent");

  mPicoDst->SetBranchStatus("*",0);         // turns OFF all branches  (speeds it up :-)
  unsigned int found;
  //mPicoDst->SetBranchStatus("EpdHit*",1,&found);   // note you need the asterisk
  //cout << "EpdHit Branch returned found= " << found << endl;
  //mPicoDst->SetBranchAddress("EpdHit",&mEpdHits);
  mPicoDst->SetBranchStatus("Event*",1,&found);
  cout << "Event Branch returned found= " << found << endl;
  mPicoDst->SetBranchAddress("Event",&mEventClonesArray);


  mPicoDst->SetBranchStatus("Track*",1,&found);
  cout << "Track Branch returned found= " << found << endl;
  mPicoDst->SetBranchAddress("Track",&mTracks);
}

//=================================================
short PicoAnalyzer::Init(char const* TPCWeightFile, char const* TPCShiftFile){
  
  // Miscellaneous One-dimensional histograms
  //mHisto1D[2] = new TH1D("dNdeta","dNdeta",100,-5.5,5.5);
  //mHisto1D[3] = new TH1D("RAWpT","RAWpT",400,0.0,4.0);

//---------------Make histograms for weighting TPC tracks----------------- 
  mNumberOfTrackTypes =mNpTbin*mNVzbin*mNEtabin*2;
  mNumberOfTrackTypes2 =mNpTbin*mNVzbin*2;
  #ifdef _OutputTPCWeightFile
  TString OutputTPCWeightFileName = mFileNameBase;
  OutputTPCWeightFileName += "_TPCWeightHists.root";
  mOutputTPCWeightFile = new TFile(OutputTPCWeightFileName,"RECREATE");

  for(int icent=0;icent<9;icent++){
    mTPCPhiWeightOutput[icent] = new TH2D(Form("TPCPhiWeightOutputCent%d",icent),Form("TPCPhiWeightOutputCent%d",icent),mNPhibin,-TMath::Pi(),TMath::Pi(),mNumberOfTrackTypes,0.5,0.5+mNumberOfTrackTypes);
    mTPCPhiAveraged[icent] = new TH2D(Form("TPCPhiAveragedCent%d",icent),Form("TPCPhiAveragedCent%d",icent),mNPhibin,-TMath::Pi(),TMath::Pi(),mNumberOfTrackTypes,0.5,0.5+mNumberOfTrackTypes);
    mTPCEtaWeightOutput[icent] = new TH2D(Form("TPCEtaWeightOutputCent%d",icent),Form("TPCEtaWeightOutputCent%d",icent),mNEtabin,mEtaMin,mEtaMax,mNumberOfTrackTypes2,0.5,0.5+mNumberOfTrackTypes2);
    mTPCEtaAveraged[icent] = new TH2D(Form("TPCEtaAveragedCent%d",icent),Form("TPCEtaAveragedCent%d",icent),mNEtabin,mEtaMin,mEtaMax,mNumberOfTrackTypes2,0.5,0.5+mNumberOfTrackTypes2);
  }
  #endif

//------------------------------------------------------------------------

//----------------Make histograms for shifting Psi_TPC--------------------
#ifdef _OutputTPCShiftFile 
  TString OutputTPCShiftFileName = mFileNameBase;
  OutputTPCShiftFileName += "_TPCShiftHists.root";
  mOutputTPCShiftFile = new TFile(OutputTPCShiftFileName,"RECREATE");
  //----------------Make histograms for shifting Psi_TPC--------------------
  for(int ivz=0;ivz<16;ivz++){
    mTPCCosShift[ivz] = new TProfile3D(Form("TPCCosShiftVz%d",ivz),Form("TPCCosShiftVz%d",ivz),9,-0.5,8.5,mFourierOrder,0.5,0.5+mFourierOrder,10,0.5,10.5);
    mTPCCosShift[ivz]->Sumw2();
    mTPCSinShift[ivz] = new TProfile3D(Form("TPCSinShiftVz%d",ivz),Form("TPCSinShiftVz%d",ivz),9,-0.5,8.5,mFourierOrder,0.5,0.5+mFourierOrder,10,0.5,10.5);
    mTPCSinShift[ivz]->Sumw2();//only consider 1st order.cent, fourier order,TPCetarange 
    
  }
#endif

//---------------Make histograms for checking the flattening of Psi_TPC----
  #ifdef _OutputTPCQAFile
  TString OutputRootFileName = mFileNameBase;
  OutputRootFileName += "_QAHists.root";
  mHistoFile = new TFile(OutputRootFileName,"RECREATE");
  mHisto1D[0] = new TH1D("Vz","Vz",100,-80,80);
  mHisto1D[1] = new TH1D("RefMult","RefMult",100,-10,600);
  for(int cent=0;cent<9;cent++){
    mPhiDisWeighted[cent] = new TH2D(Form("PhiDisWeightedCent%d",cent),Form("PhiDisWeightedCent%d",cent),mNPhibin,-TMath::Pi(),TMath::Pi(),mNumberOfTrackTypes,0.5,0.5+mNumberOfTrackTypes);
    mEtaDisWeighted[cent] = new TH2D(Form("EtaDisWeightedCent%d",cent),Form("EtaDisWeightedCent%d",cent),mNEtabin,mEtaMin,mEtaMax,mNumberOfTrackTypes2,0.5,0.5+mNumberOfTrackTypes2);
    for(int ivz=0;ivz<16;ivz++){
      mEtaPhiDisPhiWeighted[cent][ivz] = new TH2D(Form("EtaPhiDisPhiWeightedCent%dVz%d",cent,ivz),Form("EtaPhiDisPhiWeightedCent%dVz%d",cent,ivz),mNPhibin,-TMath::Pi(),TMath::Pi(),mNEtabin,mEtaMin,mEtaMax);
      mEtaPhiDisPhiEtaWeighted[cent][ivz] = new TH2D(Form("EtaPhiDisPhiEtaWeightedCent%dVz%d",cent,ivz),Form("EtaPhiDisPhiEtaWeightedCent%dVz%d",cent,ivz),mNPhibin,-TMath::Pi(),TMath::Pi(),mNEtabin,mEtaMin,mEtaMax);
      mEtaPhiDisRaw[cent][ivz] = new TH2D(Form("EtaPhiDisRawCent%dVz%d",cent,ivz),Form("EtaPhiDisRawCent%dVz%d",cent,ivz),mNPhibin,-TMath::Pi(),TMath::Pi(),mNEtabin,mEtaMin,mEtaMax);
      for(int ieta=0;ieta<10;ieta++){ 
        mTPCPsiDisWeighted[cent][ivz][ieta] = new TH1D(Form("TPC%dPsi1DisWeightedCent%dVz%d",ieta,cent,ivz),Form("TPC%dPsi1DisWeightedCent%dVz%d",ieta,cent,ivz),100,-TMath::Pi(),TMath::Pi());
        mTPCPsiDisShifted[cent][ivz][ieta] = new TH1D(Form("TPC%dPsi1DisShiftedCent%dVz%d",ieta,cent,ivz),Form("TPC%dPsi1DisShifteddCent%dVz%d",ieta,cent,ivz),100,-TMath::Pi(),TMath::Pi());
      }
    }
  }
  #endif

  #ifndef _OutputTPCWeightFile
//---------------Open histograms for getting TPCPhiWeight-----------------
  mTPCWeightFile = new TFile(TPCWeightFile,"READ");

  if(mTPCWeightFile->IsZombie()){
    std::cout<<"TPCWeightFile doesn't exist, didn't apply any phi or eta weight :("<<std::endl;
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
  #endif

  #ifdef _OutputTPCQAFile
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
  #endif
  cout << "Init done\n";
  return 0;
}

//==============================================
// some things might need resetting when there is a new run
void PicoAnalyzer::NewRun(int runId){
  //mRunCollisionSystem = WhichSystem(runId); //tend to cause error.
  mRunId = runId;
  //mRunEt += 1; 
}

//=================================================
short PicoAnalyzer::Make(int iEvent){
  //----------------- get data --------------
  mPicoDst->GetEntry(iEvent);
  StPicoEvent* event = (StPicoEvent*)((*mEventClonesArray)[0]);

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
  mHisto1D[1]->Fill(mRefMult);
  if (CentId<0) return 0;            // 80-100% - very peripheral
  mHisto1D[0]->Fill(PV.Z());

  int VzBin;
  double VzArr[17]={-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40};
  for(int i=0;i<17;i++){
    if(PV.Z()>VzArr[i]&&PV.Z()<=VzArr[i+1]){
      VzBin=i;
      break;
    }
  }
  //------------Prepare Qvectors for TPC EP (using the whole TPC)--------------
  #ifndef _OutputTPCWeightFile
  TH1D* PCosPhi;
  TH1D* PSinPhi;
  PCosPhi = new TH1D(Form("PCosPhi"),Form("PCosPhi"),10,0.5,10.5);//10 TPC eta range
  PSinPhi = new TH1D(Form("PSinPhi"),Form("PSinPhi"),10,0.5,10.5);
  
  float TPCEtaBins[21]={-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
  #endif
  //------------Begin loop over TPC tracks--------------------------
  for(int itrk=0; itrk<mTracks->GetEntries(); itrk++){
    StPicoTrack* track = (StPicoTrack*)((*mTracks)[itrk]);
    if ((track->nHitsFit()<mNhitsfit)||(!track->isPrimary())) continue; // I should check DCA too, but am confused how
    double nHitsFitRatio = track->nHitsFit()*1.0/track->nHitsMax();
    //cout<<nHitsFitRatio<<endl;
    //if (nHitsFitRatio<mNhitsfitratio) continue;
    TVector3 pMom = track->pMom();//track->gMom() if I want to look at the global tracks.
    double Tphi = pMom.Phi();
    double Teta = pMom.Eta();
    double TPt = pMom.Pt();//@pPt correct?
    double dca = track->gDCA(PV).Mag();
    if (dca>mDCAcut) continue;
    if(fabs(Teta)>=mEtaMax) continue;//for making the weighting histos, use a wider eta cut:|eta|<1.2
    if(TPt<=mpTMin||TPt>=mpTMax) continue;
    
    int TtrId=FindTrackId(track->charge(),PV.Z(),Teta,TPt);
    int TtrId2=FindTrackId2(track->charge(),PV.Z(),TPt);
    //-------TPC track weighting histos---------
    if(TtrId<=0||TtrId>mNumberOfTrackTypes){
      cout<<"Ah-oh, invalid TrackId :("<<endl;
    }
    if(TtrId2<=0||TtrId2>mNumberOfTrackTypes2){
      cout<<"Ah-oh, invalid TrackId2 :("<<endl;
    }
    
    #ifdef _OutputTPCWeightFile
    mTPCPhiWeightOutput[CentId]->Fill(Tphi,TtrId);//Fill the histogram for phi-weighting the TPC tracks
    for(int i=0;i<mNPhibin;i++){
      mTPCPhiAveraged[CentId]->Fill(-pi+2*pi*(i+0.5)/mNPhibin,TtrId,1.0/mNPhibin);
    }
    mTPCEtaWeightOutput[CentId]->Fill(Teta,TtrId2);//Fill the histogram for phi-weighting the TPC tracks
    mTPCEtaAveraged[CentId]->Fill(mTPCEtaAveraged[CentId]->GetXaxis()->GetBinCenter(mTPCEtaAveraged[CentId]->GetXaxis()->FindBin(Teta)),TtrId2,0.5);
    mTPCEtaAveraged[CentId]->Fill(mTPCEtaAveraged[CentId]->GetXaxis()->GetBinCenter(mTPCEtaAveraged[CentId]->GetXaxis()->FindBin(-Teta)),TtrId2,0.5);
    #endif
     
    #ifndef _OutputTPCWeightFile
    double TrackPhiWeight=1.0;
    if(mTPCPhiWeightInput[CentId]!=0) TrackPhiWeight=1.0/mTPCPhiWeightInput[CentId]->GetBinContent(mTPCPhiWeightInput[CentId]->FindBin(Tphi,TtrId));
    if(!TMath::Finite(TrackPhiWeight)) continue;
    if(TrackPhiWeight>10.0) continue;
    double TrackEtaWeight=1.0;
    if(mTPCEtaWeightInput[CentId]!=0) TrackEtaWeight=1.0/mTPCEtaWeightInput[CentId]->GetBinContent(mTPCEtaWeightInput[CentId]->FindBin(Teta,TtrId2));
    if(!TMath::Finite(TrackEtaWeight)) continue; 
    double TrackWeight=TrackPhiWeight*TrackEtaWeight;  
    if(!TMath::Finite(TrackWeight)) continue;

    #ifdef _OutputTPCQAFile
    mPhiDisWeighted[CentId]->Fill(Tphi,TtrId,TrackPhiWeight);
    mEtaDisWeighted[CentId]->Fill(Teta,TtrId2,TrackEtaWeight);
    mEtaPhiDisRaw[CentId][VzBin]->Fill(Tphi,Teta);
    mEtaPhiDisPhiWeighted[CentId][VzBin]->Fill(Tphi,Teta,TrackPhiWeight);
    mEtaPhiDisPhiEtaWeighted[CentId][VzBin]->Fill(Tphi,Teta,TrackWeight);
    #endif

    if(TMath::Abs(Teta)>=mEtaMaxAly) continue;//For the following EP measurement |eta|<1.0
    double TPCWeight[_PsiOrderMax];
    TPCWeight[0] = TrackWeight*(-1.0*Teta);//weight the TPC tracks with (-eta) for Psi1
    //TPCWeight[1] = TrackPhiWeight*TPt;//weight the TPC tracks with pT for Psi2

    //10 TPC eta range  
    for(int ieta=0;ieta<10;ieta++){
      if(fabs(Teta)<(1.0-0.1*(double)ieta)){
        PCosPhi->Fill(ieta+1,TPCWeight[0]*cos(Tphi));
        PSinPhi->Fill(ieta+1,TPCWeight[0]*sin(Tphi));
      }
    }
    #endif
    
  } 
//-------------End loop over TPC tracks-----------------------

//-------------Now let's play with the EP---------------------
  #ifndef _OutputTPCWeightFile
  double TPCPsi[10];//0_[-1.0,1.0],1_[-0.9,0.9],2_[-0.8,0.8]...9_[-0.1,0.1]
  for(int i=1; i<=10;i++){
    TPCPsi[i-1]=-999;
    if(fabs(PSinPhi->GetBinContent(i))>1e-6&&fabs(PCosPhi->GetBinContent(i))>1e-6){
      TPCPsi[i-1] = atan2(PSinPhi->GetBinContent(i),PCosPhi->GetBinContent(i));
    }
    //atan2 returns angle between -pi and pi, remember to divide by i!!!! 06/16/2020
  }
  
  #ifdef _OutputTPCShiftFile
  //--------Fill the histos for Psi-shifting the TPC_EP---------
  for(int k=1;k<=mFourierOrder;k++){
    for(int i=1;i<=10;i++){
      double tmp=(double)k;
      if(TPCPsi[i-1]>-100){
        mTPCCosShift[VzBin]->Fill(CentId,k,i,cos(tmp*TPCPsi[i-1]));
        mTPCSinShift[VzBin]->Fill(CentId,k,i,sin(tmp*TPCPsi[i-1]));
      }
    }
  }
  #endif

  #ifdef _OutputTPCQAFile
  //-------------Shift the TPC EP if the shifting histos exist------------
  double TPCPsiShifted[10];
  for(int i=0;i<10;i++){
    TPCPsiShifted[i]=TPCPsi[i];
  }
  if(mTPCCosShift_Input[VzBin]!=0&&mTPCSinShift_Input[VzBin]!=0){
    double deltapsi[10]={0.0};
    for(int i=1;i<=10;i++){
      for(int k=1;k<=mFourierOrder;k++){
        double tmp=(double)(1.0*k);
        deltapsi[i-1]+=2.0*(-mTPCSinShift_Input[VzBin]->GetBinContent(CentId+1,k,i)*cos(tmp*TPCPsi[i-1])+mTPCCosShift_Input[VzBin]->GetBinContent(CentId+1,k,i)*sin(tmp*TPCPsi[i-1]))/tmp;
      }
      if(TPCPsiShifted[i-1]>-100){
        TPCPsiShifted[i-1]+=deltapsi[i-1];
        if(TPCPsiShifted[i-1]<-pi) TPCPsiShifted[i-1]+=2.0*pi;
        if(TPCPsiShifted[i-1]>pi) TPCPsiShifted[i-1]-=2.0*pi;
      }
    }  
  }


  //----Fill the histograms for checking the flattening of the TPC EPs-----
  for(int i=0;i<10;i++){
    if(TPCPsi[i]>-100&&TPCPsiShifted[i]>-100){
      mTPCPsiDisWeighted[CentId][VzBin][i]->Fill(TPCPsi[i]);
      mTPCPsiDisShifted[CentId][VzBin][i]->Fill(TPCPsiShifted[i]);
    }
  }
  #endif

  delete PCosPhi;
  delete PSinPhi;
  #endif
  return 0;
}
//=================================================
short PicoAnalyzer::Finish(){

  #ifndef _OutputTPCWeightFile
  mTPCWeightFile->Close();
  #ifdef _OutputTPCQAFile
  mTPCShiftFile->Close();
  #endif
  #endif

  #ifdef _OutputTPCWeightFile
  mOutputTPCWeightFile->Write();
  mOutputTPCWeightFile->Close();
  #endif
  
  #ifdef _OutputTPCShiftFile
  mOutputTPCShiftFile->Write();
  mOutputTPCShiftFile->Close();
  #endif

  #ifdef _OutputTPCQAFile
  mHistoFile->Write();
  mHistoFile->Close();
  #endif

  cout << "Finish!!\n\n";

  return 0;
}

//TrackId starts at 1!
int PicoAnalyzer::FindTrackId(int Trkch,double TrkVz,double TrkEta,double TrkpT){
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

int PicoAnalyzer::FindTrackId2(int Trkch,double TrkVz,double TrkpT){
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


