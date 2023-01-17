#ifndef PicoAnalyzer__
#define PicoAnalyzer__

#include "TObject.h"

class TChain;
class TClonesArray;
class StPicoEvent;
class TH1D;
class TH2D;
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
class StRefMultCorr;
#include "StRefMultCorr/CentralityMaker.h"
#include "TString.h"
#include "TMath.h"

#define _PsiOrderMax 2 //Maximum order of EP to worry about 

class PicoAnalyzer : public TObject {
 public:
  PicoAnalyzer(TString FileNameBase="MikesStuff");
  ~PicoAnalyzer();

  void SetPicoDst(TChain*);
  short Init(char const* TPCweightFileName="TPCWeightFile.root",char const*  TPCShiftFileName="TPCShiftFile.root");
  short Make(int iEvent);
  short Finish();

 private:

  TString mFileNameBase;

  // parameters relevant to my analysis
  //double mNmipQtB;  // ADC value of MIP peak on rings 6-16 (read out thru QT32Bs)
  //double mNmipQtC;  // ADC value of MIP peak on rings 1-5  (read out thru QT32Cs)
  //double mnMipThreshold;  // low-signal threshold, to cut out noise basically.

  // the data objects
  TChain*   mPicoDst;
  //TClonesArray* mEpdHits;
  //TClonesArray* mBbcHits;
  TClonesArray* mTracks;
  TClonesArray* mEventClonesArray;  // kind of hilarious that the StPicoEvent is stored as a one-element TClonesArray :-)

  int mRunId;                         // when this changes, refresh some information.
  int mRunEt;                    //Run entry
  short mRunCollisionSystem;

  // internal methods
  int FindTrackId(int Trkch,double TrkVz,double Trketa,double TrkpT );//this method only works when bin numbers are even and non-zero!!!!
  int FindTrackId2(int Trkch,double TrkVz,double TrkpT );//this method only works when bin numbers are even and non-zero!!!!
  void ReadInSystems();    // reads a text file that idenfies the collision system for every run
  short WhichSystem(int runId);
  void NewRun(int runId);    // invoked when Make() detects that a new run has been loaded
  //void FillPhiWeightHistos(StPicoEpdHit* epdHit, double weight);     // fills the histograms used for (a later job's) phi weighting
  bool Runlist(int runId);
  // https://drupal.star.bnl.gov/STAR/blog/lisa/optimizing-ep1-resolution-au27au-ring-dependent-weights-eta-dependent-weights
  double v1Weight(int CentId, double eta);

  // useful objects kept by PicoAnalyzer
  //StEpdGeom* mEpdGeom;
  //StBbcGeom* mBbcGeom;
  //StEpdEpFinder* mEpFinder;
  StRefMultCorr* mRefMultCorr;
  TRandom3* mRan;//seems like RCF only likes TRandom3
  char mCollidingSystem[365][500];    // index1=day of year;  index2=run of day
  
  //static const int mTPCphibin = 80;
  //double mTPCphibinwidth = TMath::Pi()/40;
  

  // ---------- Now, my histograms, ntuples, TFiles, etc.  All stuff particular to my analysis
  // 1D histograms
  TH1D* mHisto1D[8];            // miscellaneous 1D histograms
  
  // 2D histograms
  //TH2D* mHisto2D[40];           // miscellaneous 2D histograms
  TH2D* mEPDPsiEvsW[_PsiOrderMax][3];//row, weighted and shifted

  TH2D* mTPCPhiWeightOutput[9];    // this is used for "phi weighting" TPC tracks 
  TH2D* mTPCPhiAveraged[9];        
  TH2D* mTPCPhiWeightInput[9];     // "phi weighting" correction factors that were calculated and saved in a PREVIOUS run
  TH2D* mTPCPhiAveragedInput[9];
  
  TH2D* mTPCEtaWeightOutput[9];    // this is used for "eta weighting" TPC tracks 
  TH2D* mTPCEtaAveraged[9];        
  TH2D* mTPCEtaWeightInput[9];     // "eta weighting" correction factors that were calculated and saved in a PREVIOUS run
  TH2D* mTPCEtaAveragedInput[9];

double mpTMin;
double mpTMax;
double mEtaMin;
double mEtaMax;
double mVzMin;
double mVzMax;
int mNpTbin;
int mNVzbin;
int mNEtabin;
int mNPhibin;
int mNumberOfTrackTypes;
int mNumberOfTrackTypes2;
double mVtxR;
double mDiffVzVPD;
int mNhitsfit;
double mNhitsfitratio;
double mDCAcut;
int mFourierOrder;
double mEtaMaxAly;

//histos for checking the phi-weighted eta-phi distribution for the TPC tracks 
  TH2D* mEtaPhiDisPhiWeighted[9][16];//phi-weighted eta phi distribution for 9 centralitites.
  TH2D* mEtaPhiDisPhiEtaWeighted[9][16];//phi-weighted eta phi distribution for 9 centralitites.
  TH2D* mEtaPhiDisRaw[9][16];
  TH2D* mPhiDisWeighted[9]; 
  TH2D* mEtaDisWeighted[9]; 
//histos for checking the pT distribution 
  //TH1F* mPtRaw[9][10];//9 centralities and 10 eta windows
  //TH1D* mpTPhiWeighted[9][10];
  TH1F* mPtRaw;

//histos for checking the flattening of the TPC Psi
  TH1D* mTPCPsiDisWeighted[9][16][10];//Entry is the centrality
  TH1D* mTPCPsiDisShifted[9][16][10];
  
  //TH1D* mEPDFullPsiWeighted[9][_PsiOrderMax];
  //TH1D* mEPDFullPsiShifted[9][_PsiOrderMax];

//TProfile for calculating the resolutions 
  //TProfile* mResolution[9][16][10];

  TProfile3D* mTPCCosShift[16];
  TProfile3D* mTPCSinShift[16];

  TProfile3D* mTPCCosShift_Input[16];
  TProfile3D* mTPCSinShift_Input[16];


  //TProfile2D* mProfile2D[40];    // miscellaneous 2D profiles


  // TFiles (to store histograms and data)
  TFile* mHistoFile;
  TFile* mTPCWeightFile;
  TFile* mTPCShiftFile;//input file for psi shifting in "run2"
  TFile* mOutputTPCWeightFile;
  TFile* mOutputTPCShiftFile;

  static const int mEPTPCMaxTerm = 6; 
  ClassDef(PicoAnalyzer, 1)                     //  Macro for CINT compatability

};


#endif
