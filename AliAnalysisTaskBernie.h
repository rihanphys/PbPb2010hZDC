/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */
/* $Id: $ */

#ifndef AliAnalysisTaskBernie_H
#define AliAnalysisTaskBernie_H

/////////////////////////////////////////////////
// AliAnalysisTaskBernie:
// analysis task for Scalar Product method
// Author: Naomi van der Kolk (kolk@nikhef.nl)
/////////////////////////////////////////////////

class AliFlowEventSimple;
//class AliFlowAnalysisIDCSP;
class TList;

#include "TString.h"
#include "AliAnalysisTaskSE.h"
#include "TNtuple.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile2D.h"
#include "TProfile.h"

//===============================================================

class AliAnalysisTaskBernie : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskBernie();
  AliAnalysisTaskBernie(const char *name, Bool_t usePhiWeights);
  virtual ~AliAnalysisTaskBernie();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void   SetUsePhiWeights(Bool_t const aPhiW) {this->fUsePhiWeights = aPhiW;}
  Bool_t GetUsePhiWeights() const             {return this->fUsePhiWeights;}

  void     SetRelDiffMsub(Double_t diff) { this->fRelDiffMsub = diff; }
  Double_t GetRelDiffMsub() const        { return this->fRelDiffMsub; }

  void SetApplyCorrectionForNUA(Bool_t const applyCorrectionForNUA) {this->fApplyCorrectionForNUA = applyCorrectionForNUA;};
  Bool_t GetApplyCorrectionForNUA() const {return this->fApplyCorrectionForNUA;};

  void SetHarmonic(Int_t const harmonic) {this->fHarmonic = harmonic;}
  Int_t GetHarmonic() const {return this->fHarmonic;};

  void SetBehaveAsEP() { fNormalizationType = 0; }
  void SetTotalQvector(const char *tqv) { *this->fTotalQvector = tqv; }
  void SetBookOnlyBasicCCH(Bool_t const aMB) { this->fMinimalBook = aMB; }

  //rihan:
  void SetRunFlag(Int_t const runnum) { this->frunflag = runnum; }
  void SetZDCESEList(TList* const kList) {this->fZDCESEList = kList;};

  
private:

  AliAnalysisTaskBernie(const AliAnalysisTaskBernie& aAnalysisTask);
  AliAnalysisTaskBernie& operator=(const AliAnalysisTaskBernie& aAnalysisTask);

  AliFlowEventSimple* fEvent;         //input event

//AliFlowAnalysisIDCSP** fSP;         // analysis object
//AliFlowAnalysisIDCSP*  fSP;         // analysis object only one.
  Int_t              fievent;
  TList         *fListHistos;         // collection of output
  TList        *fListWeights;         // list with weights
  Bool_t        fMinimalBook;         // flag to turn off QA and minimize FlowCommonHist
  Bool_t      fUsePhiWeights;         // use phi weights
  Int_t            fHarmonic;         // harmonic
  Int_t   fNormalizationType;         // 0: EP mode || 1: SP mode (default)
  Int_t           fNCentBins;         // Number of Cent bins
  Double_t      fRelDiffMsub;         // the relative difference the two subevent multiplicities can have
  TString     *fTotalQvector;         // total Q-vector is: "QaQb" (means Qa+Qb), "Qa"  or "Qb"
  Bool_t fApplyCorrectionForNUA;      // apply automatic correction for non-uniform acceptance
  
  //rihan:
  Int_t             frunflag; 
  Int_t            oldRunNum; 
  Int_t          runNums[90];  

  Float_t           VxCut[2];  //
  Float_t           VyCut[2];  //
  Float_t          maxVx_low;
  Float_t         maxVx_high;
  Float_t          maxVy_low; 
  Float_t         maxVy_high;
  Int_t            checkOnce;  //
  TList         *fZDCESEList;  

  TProfile2D  *fHist_znCx_V0_VxVy[89];   //!
  TProfile2D  *fHist_znCy_V0_VxVy[89];   //!
  TProfile2D  *fHist_znAx_V0_VxVy[89];   //!
  TProfile2D  *fHist_znAy_V0_VxVy[89];   //!

  TProfile    *fHist_Vx_vs_runnumber;    //!
  TProfile    *fHist_Vy_vs_runnumber;    //!

  TProfile2D        *fHist_QnxRecent;    //!
  TProfile2D        *fHist_QnyRecent;    //!

  TH1F              *fHist_vBincount;    //!
  TH1F              *fHist_Psi1_zdnA;    //!
  TH1F              *fHist_Psi1_zdnC;    //!
  TH1F            *fHist_EventperRun;    //!
  
  TH1F            *fHist_EventPlane2;    //!
  TH1F         *fHist_Vx_ArrayFinder;    //!
  TH1F         *fHist_Vy_ArrayFinder;    //!
  TH1F         *fHist_Vz_ArrayFinder;    //!
  TH1F              *fHist_Vertex_Vz;    //!
  TH1F            *fHist_Event_count;    //!
  TH1F            *fHist_Centr_count;    //!
  TH2F               *fHistQA_etaphi;    //!
  
  TH2F              *fHist_Vertex_XY;    //!
  TH2F              *fHist_ZDCn_A_XY;    //!
  TH2F              *fHist_ZDCn_C_XY;    //!
  TProfile       *fHist_Vx_vs_runnum;    //!
  TProfile       *fHist_Vy_vs_runnum;    //!
  TProfile       *fHist_Vz_vs_runnum;    //!
  TProfile   *fHist_tracks_vs_runnum;    //!

  
  ClassDef(AliAnalysisTaskBernie, 2); // example of analysis
};

//==================================================================

#endif
