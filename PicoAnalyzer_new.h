#ifndef PicoAnalyzer_new__
#define PicoAnalyzer_new__


#include "TObject.h"

class TChain;
class TClonesArray;
class StPicoEvent;
class TH1D;
class TH2D;
class TH3D;
class TH3F;
class TH2F;
class TH1F;
class TProfile;
class TProfile2D;
class TProfile3D;
class TNtuple;
class TFile;
class StEpdGeom;
class StBbcGeom;
class TRandom3;
class StPicoEpdHit; 
class StEpdEpFinder;
class TLorentzVector;
class StRefMultCorr;
#include "StRefMultCorr/CentralityMaker.h"
#include "TString.h"
#include "TMath.h"

#define _PsiOrderMax 1 //Maximum order of EP to worry about 

class PicoAnalyzer_new : public TObject {
 public:
  PicoAnalyzer_new(TString FileNameBase="MikesStuff");
  ~PicoAnalyzer_new();

  void SetPicoDst(TChain*);
  //short Init(char const* TPCweightFileName="TPCWeightFile.root",char const*  TPCShiftFileName="TPCShiftFile.root",char const* EPDPhiWeightFile="EPDcorrection.root");
  short Init(char const* TPCweightFileName="TPCWeightFile.root",char const*  TPCShiftFileName="TPCShiftFile.root");
  short Make(int iEvent);
  short Finish();

 private:

  TString mFileNameBase;
  
  // the data objects
  TChain*   mPicoDst;
  TClonesArray* mEpdHits;
  TClonesArray* mBbcHits;
  TClonesArray* mTracks;
  TClonesArray* mEventClonesArray;  // kind of hilarious that the StPicoEvent is stored as a one-element TClonesArray :-)
  TClonesArray* mTraits;

  int mRunId;                         // when this changes, refresh some information.
  int mRunEt;                    //Run entry
  short mRunCollisionSystem;


  // internal methods
  int FindTrackId(int Trkch,double TrkVz,double Trketa,double TrkpT );//this method only works when bin numbers are even and non-zero!!!!
  int FindTrackId2(int Trkch,double TrkVz,double TrkpT );//this method only works when bin numbers are even and non-zero!!!!
  int FindCent(int RefMult);   // utility class just giving centrality bin.  Copied directly from Isaac 1 May 2018
  double GetBbcPmtPhi(short PmtId);
  void ReadInSystems();    // reads a text file that idenfies the collision system for every run
  short WhichSystem(int runId);
  void NewRun(int runId);    // invoked when Make() detects that a new run has been loaded
  //void FillPhiWeightHistos(StPicoEpdHit* epdHit, double weight);     // fills the histograms used for (a later job's) phi weighting
  bool Runlist(int runId);
  // https://drupal.star.bnl.gov/STAR/blog/lisa/optimizing-ep1-resolution-au27au-ring-dependent-weights-eta-dependent-weights
  double v1Weight(int CentId, double eta);

  // useful objects kept by PicoAnalyzer_new
  StEpdGeom* mEpdGeom;
  StBbcGeom* mBbcGeom;
  StEpdEpFinder* mEpFinder;
  StRefMultCorr* mRefMultCorr;
  TRandom3* mRan;//seems like RCF only likes TRandom3
  char mCollidingSystem[365][500];    // index1=day of year;  index2=run of day
  
  static const int mTPCphibin = 80;
  double mTPCphibinwidth = TMath::Pi()/40;
  
  //TH2D* mTPCPhiWeightOutput;    // this is used for "phi weighting" TPC tracks 
  //TH2D* mTPCPhiAveraged;        // this is just used for the normalization of the above. - not ever saved.  Only internal use
  TH2D* mTPCPhiWeightInput[9];     // "phi weighting" correction factors that were calculated and saved in a PREVIOUS run
  TH2D* mTPCPhiAveragedInput[9];
  TH2D* mTPCEtaWeightInput[9];     // "eta weighting" correction factors that were calculated and saved in a PREVIOUS run
  TH2D* mTPCEtaAveragedInput[9];

  double mpTMin;
  double mpTMax;
  double mEtaMin;//for phi-weighting, used in FindTrackId() and when MakeWeight
  double mEtaMax;////for phi-weighting
  int mNPhibin;
  double mVzMin;
  double mVzMax;
  double mEPDMax;
  double mEPDthresh;
  int mNpTbin;
  int mNVzbin;
  int mNEtabin;
  int mNumberOfTrackTypes;
  int mNumberOfTrackTypes2;
  double mVtxR;
  double mDiffVzVPD;
  int mNhitsfit;
  double mNhitsfitratio;
  double mDCAcut;
  int mFourierOrder;
  int mNTPCSubEvents;
  int mNEPDSubEvents;
  double mpTbound;//For looking at low/high pT tracks
  double mPionSigma;//1/beta sigma 
  double mKaonSigma;//1/beta sigma
  double mProtonSigma;//1/beta sigma 
  double mEtaMaxv1;//for v1 analysis
  double mEtaMinv1;//for v1 analysis 

//histos for checking the phi-weighted eta-phi distribution for the TPC tracks 
  TH2D* mEtaPhiDisPhiWeighted[9];//phi-weighted eta phi distribution for 9 centralitites.
  TH2D* mEtaPhiDisRaw[9];
//histos for checking the pT distribution 
  TH1F* mPtRaw[9];//9 centralities 
  //TH1D* mpTPhiWeighted[9][10];
  //TH1F* mPtRaw;
  TH1D* mVz[9];//Vz dis for nine centralites
  TH1D* mdNdeta[9];//dNdeta for nine centralities 

//histos for checking the flattening of the TPC Psi
  TH1F* mTPCPsiDisWeighted[9][16][_PsiOrderMax];//9 cent; 16 vz bins;
  TH1F* mTPCPsiDisShifted[9][16][_PsiOrderMax];
//histos for checking the flattening of the Psi_EODfull
  TH1F* mEPDFullPsiWeighted[9][16][_PsiOrderMax];
  TH1F* mEPDFullPsiShifted[9][16][_PsiOrderMax];
  TH2F* mEPDPsiEvsWWeighted[9][16][_PsiOrderMax];
  TH2F* mEPDPsiEvsWShifted[9][16][_PsiOrderMax];

//TProfile for calculating the resolutions.
  TProfile* mResolution[9][16][_PsiOrderMax];//Indicies are cent. 


//TH3D for measuring dNdphi
  //TH3D* mThreeD[9][2];//centralities, ew
  TH3F* mThreeDTile[9][16][2][24];//cent, vz bin, ew, tile 

  TProfile3D* mTPCCosShift[16];
  TProfile3D* mTPCSinShift[16];

  TProfile3D* mTPCCosShift_Input[16];
  TProfile3D* mTPCSinShift_Input[16];

  TProfile2D* mAveEta[9][2]; //nine centralities, EW(ring#, ivz,eta)


  // TFiles (to store histograms and data)
  TFile* mHistoFile;
  TFile* mTPCWeightFile;//input 
  TFile* mTPCShiftFile;//input
  TFile* mOutputTPCShiftFile;

  static const int mEPTPCMaxTerm = 6; 
  ClassDef(PicoAnalyzer_new, 1)                     //  Macro for CINT compatability

};


#endif