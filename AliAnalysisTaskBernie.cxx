/*************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

///////////////////////////////////////////////
// AliAnalysisTaskBernie:
// rihan: there was a task for Scalar product. Bernie edited it.
// Author: Naomi van der Kolk (kolk@nikhef.nl)
///////////////////////////////////////////////


#include "Riostream.h" //needed as include

class TFile;
class TList;
class AliAnalysisTaskSE;

#include "TProfile.h"  //needed as include
#include "TNtuple.h"
#include "TProfile2D.h"
#include "TH1.h"
#include "TH2.h"
#include "TApplication.h"
#include "TTree.h"

#include "AliAnalysisManager.h"
#include "AliFlowEventSimple.h"
#include "AliAnalysisTaskBernie.h"
#include "AliFlowCommonHist.h"
#include "AliFlowCommonHistResults.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliESDInputHandler.h"
#include "AliESDZDC.h"
#include "AliMultiplicity.h"
#include "AliAnalysisUtils.h"
#include "AliAODHandler.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODVertex.h"
#include "AliAODVZERO.h"
#include "AliAODZDC.h"
#include "AliAODMCHeader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliHeader.h"
#include "AliVParticle.h"
#include "AliStack.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliGenEventHeader.h"
#include "AliPhysicsSelectionTask.h"
#include "AliPhysicsSelection.h"
#include "AliBackgroundSelection.h"
#include "AliTriggerAnalysis.h"
#include "AliCentrality.h"


#include "AliLog.h"

using std::endl;
using std::cout;
ClassImp(AliAnalysisTaskBernie)

//________________________________________________________________________
AliAnalysisTaskBernie::AliAnalysisTaskBernie(const char *name, Bool_t usePhiWeights) :
AliAnalysisTaskSE(name),
fEvent(NULL),
//fSP(NULL),
fListHistos(NULL),
fMinimalBook(kFALSE),
fUsePhiWeights(usePhiWeights),
fListWeights(NULL),
fRelDiffMsub(1.0),
fApplyCorrectionForNUA(kFALSE),
fHarmonic(2),
fNormalizationType(1),
fNCentBins(9),
fTotalQvector(NULL),
fievent(0),
frunflag(0),
fHist_Vertex_XY(NULL),
fHist_Vx_vs_runnum(NULL),
fHist_Vy_vs_runnum(NULL),
fHist_Vz_vs_runnum(NULL),
fHist_tracks_vs_runnum(NULL),
fHist_Vx_ArrayFinder(NULL),
fHist_Vy_ArrayFinder(NULL),
fHist_Vz_ArrayFinder(NULL),
fHist_Vertex_Vz(NULL),
fHist_Event_count(NULL),
fHist_Centr_count(NULL),
fHist_ZDCn_A_XY(NULL),
fHist_ZDCn_C_XY(NULL),
fHist_EventPlane2(NULL),
fHist_QnxRecent(NULL),
fHist_QnyRecent(NULL),  
fHist_EventperRun(NULL),
fHist_Psi1_zdnA(NULL),
fHist_Psi1_zdnC(NULL),
fHist_vBincount(NULL),
fHistQA_etaphi(NULL),
fZDCESEList(NULL),
fHist_Vx_vs_runnumber(NULL),
fHist_Vy_vs_runnumber(NULL),
oldRunNum(0),
maxVx_low(0),
maxVx_high(0),
maxVy_low(0), 
maxVy_high(0),
checkOnce(0)
{
  //AliDebug(2,"AliAnalysisTaskBernie::AliAnalysisTaskBernie(const char *name)");

   DefineInput(1, AliFlowEventSimple::Class()); //Input slot #0 works with an AliFlowEventSimple

   for(int i=0;i<frunflag;i++)
     runNums[i] = 0;
   //printf("\n ...  AliAnalysisTaskBernie() constructor called 3...  \n");  

  
for(int i=0;i<89;i++){
       fHist_znCx_V0_VxVy[i] = NULL;
       fHist_znCy_V0_VxVy[i] = NULL;
       fHist_znAx_V0_VxVy[i] = NULL;
       fHist_znAy_V0_VxVy[i] = NULL;
 }

     for(int i=0;i<2;i++)
    {
      VxCut[i] = 0;
      VyCut[i] = 0;
    }
   
           DefineOutput(1,TList::Class());
	   
	  //Input slot #2 is needed for the weights input file
	  if(usePhiWeights) {
	    DefineInput(2, TList::Class());                 // I need this for pass 2, when I have weights
	  }

	//Total Q-vector is: "QaQb" (means Qa+Qb), "Qa"  or "Qb"
	fTotalQvector = new TString("QaQb");

	printf("\n ...  AliAnalysisTaskBernie() constructor called ...  \n");
}

//________________________________________________________________________
AliAnalysisTaskBernie::AliAnalysisTaskBernie() :
AliAnalysisTaskSE(),
fEvent(NULL),
//fSP(NULL),
fListHistos(NULL),
fMinimalBook(kFALSE),
fUsePhiWeights(kFALSE),
fListWeights(NULL),
fRelDiffMsub(1.0),
fApplyCorrectionForNUA(kFALSE),
fHarmonic(0),
fNormalizationType(1),
fNCentBins(9),
fTotalQvector(NULL),
fievent(0),
frunflag(0),
fHist_Vertex_XY(NULL),
fHist_Vx_vs_runnum(NULL),
fHist_Vy_vs_runnum(NULL),
fHist_Vz_vs_runnum(NULL),
fHist_tracks_vs_runnum(NULL),
fHist_Vx_ArrayFinder(NULL),
fHist_Vy_ArrayFinder(NULL),
fHist_Vz_ArrayFinder(NULL),
fHist_Vertex_Vz(NULL),
fHist_Event_count(NULL),
fHist_Centr_count(NULL),
fHist_ZDCn_A_XY(NULL),
fHist_ZDCn_C_XY(NULL),
fHist_EventPlane2(NULL),
fHist_QnxRecent(NULL),
fHist_QnyRecent(NULL),
fHist_EventperRun(NULL),
fHist_Psi1_zdnA(NULL),
fHist_Psi1_zdnC(NULL),
fHist_vBincount(NULL),
fHistQA_etaphi(NULL),
fZDCESEList(NULL),
fHist_Vx_vs_runnumber(NULL),
fHist_Vy_vs_runnumber(NULL),
oldRunNum(0),
maxVx_low(0),
maxVx_high(0),
maxVy_low(0), 
maxVy_high(0),
checkOnce(0)
{

  
for(int i=0;i<89;i++){
       fHist_znCx_V0_VxVy[i] = NULL;
       fHist_znCy_V0_VxVy[i] = NULL;
       fHist_znAx_V0_VxVy[i] = NULL;
       fHist_znAy_V0_VxVy[i] = NULL;
 }

  for(int i=0;i<2;i++)
    {
      VxCut[i] = 0;
      VyCut[i] = 0;
    }
  
	// Constructor
	// AliDebug(2,"AliAnalysisTaskBernie::AliAnalysisTaskBernie()");
}

//________________________________________________________________________
AliAnalysisTaskBernie::~AliAnalysisTaskBernie()
{
    //delete [] fSP;          //is this not deleting all fSP[i] ?
     printf("\n ... ~ AliAnalysisTaskBernie() destructor called ...  \n");
}

//________________________________________________________________________
void AliAnalysisTaskBernie::UserCreateOutputObjects()
{
 //Called at every worker node to initialize
 AliDebug(2,"AliAnalysisTaskBernie::CreateOutputObjects() is called....");

  //run:2010
  Int_t runArray[89] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137243, 137236, 137235, 137232, 137231, 137162, 137161};
  //run:2011
  //Int_t runArray[68] = {167915, 168115, 168460, 169035, 169238, 169859, 170228, 167920, 168310, 168464, 169091, 169411, 169923, 170230, 167985, 168311, 168467, 169094, 169415, 170027, 170268, 167987, 168322, 168511, 169138, 169417, 170081, 170269, 167988, 168325, 168512, 169144, 169835, 170155, 170270, 168069, 168341, 168514, 169145, 169837, 170159, 170306, 168076, 168342, 168777, 169148, 169838, 170163, 170308, 168105, 168361, 168826, 169156, 169846, 170193, 170309, 168107, 168362, 168988, 169160, 169855, 170203, 168108, 168458, 168992, 169167, 169858, 170204};

  for(int i=0;i<89;i++)
  {
    runNums[i] = runArray[i];
  }
  
  Int_t totalQ = 0;
  if(fTotalQvector->Contains("Qa")) totalQ += 1;
  if(fTotalQvector->Contains("Qb")) totalQ += 2;

  fListHistos = new TList();
  fListHistos->SetOwner(kTRUE);
  

  fHist_vBincount      = new TH1F("fHist_vBincount"," ",500,0,500);
  fHist_EventperRun    = new TH1F("fHist_EventperRun","", frunflag,0,frunflag);
  fHist_Psi1_zdnA      = new TH1F("fHist_Psi1_zdnA","",  200, 0,2.*TMath::Pi());
  fHist_Psi1_zdnC      = new TH1F("fHist_Psi1_zdnC","",  200, 0,2.*TMath::Pi());

  fHist_EventPlane2    = new TH1F("fHist_EventPlane_TPC","",100, 0,TMath::Pi());
  fHist_Vertex_Vz      = new TH1F("fHist_Vertext_Vz_dist","", 200,    -15,    15);
  fHistQA_etaphi     = new       TH2F("fHistQA_etaphi","phi vs eta",400,0,6.4,100,-1,1); 
  fHist_ZDCn_A_XY    = new       TH2F("fHist_ZDCn_A_XY","ZDC A XY Distribution",80,-2,2,80,-2,2);
  fHist_ZDCn_C_XY    = new       TH2F("fHist_ZDCn_C_XY","ZDC C XY Distribution",80,-2,2,80,-2,2);
  
  fHist_Vertex_XY    = new       TH2F("fHist_Vertex_XY","Vertex XY Distribution",100,-0.2,0.2,100,0.0,0.4);
  fHist_Vx_vs_runnum = new TProfile("fHist_Vx_vs_runnum","<Vx>_vs_runnum",frunflag,0,frunflag);
  fHist_Vy_vs_runnum = new TProfile("fHist_Vy_vs_runnum","<Vy>_vs_runnum",frunflag,0,frunflag);
  fHist_Vz_vs_runnum = new TProfile("fHist_Vz_vs_runnum","<Vz>_vs_runnum",frunflag,0,frunflag); 
  fHist_tracks_vs_runnum = new TProfile("fHist_tracks_vs_runnum","<nTracks>_vs_runnum",frunflag,0,frunflag);

  
  fHist_Event_count    = new TH1F("fHist_Event_count"," ",20,0,20);
  
  Double_t centRange[12]    = {0,5,10,20,30,40,50,60,70,80,90,100};
  fHist_Centr_count         = new TH1F("fHist_Centr_count"," ",11,centRange);

 fListHistos->Add(fHist_vBincount);
 fListHistos->Add(fHist_EventperRun);
 fListHistos->Add(fHist_Psi1_zdnA);
 fListHistos->Add(fHist_Psi1_zdnC);
 fListHistos->Add(fHist_EventPlane2);
 fListHistos->Add(fHistQA_etaphi);
 fListHistos->Add(fHist_ZDCn_A_XY);
 fListHistos->Add(fHist_ZDCn_C_XY);  
 fListHistos->Add(fHist_Event_count);
 fListHistos->Add(fHist_Centr_count);
 fListHistos->Add(fHist_Vertex_Vz);
 fListHistos->Add(fHist_Vertex_XY);
 fListHistos->Add(fHist_Vx_vs_runnum);
 fListHistos->Add(fHist_Vy_vs_runnum);
 fListHistos->Add(fHist_Vz_vs_runnum);
 fListHistos->Add(fHist_tracks_vs_runnum);

 // random comment 
 //another one
 fHist_QnxRecent = new TProfile2D(Form("fHist_Q2xRecent"),"",90,0,90,frunflag,0,frunflag);
 fHist_QnyRecent = new TProfile2D(Form("fHist_Q2yRecent"),"",90,0,90,frunflag,0,frunflag);

 //fListHistos->Add(fHist_QnxRecent);
 //fListHistos->Add(fHist_QnyRecent);
 

 //-------- define and add the recenter histograms ----------------
  const int vxBin = 25;
  const int vyBin = 25;
  const int vzBin = 25;
  //const int NbinVt =  (vzBin-1)*vyBin*vxBin + (vyBin-1)*vxBin + vxBin;
  const int NbinVt =   (vyBin-1)*vxBin + vxBin;
   
  VxCut[0] =  -0.035;
  VxCut[1] =   0.020;
  VyCut[0] =   0.144;
  VyCut[1] =   0.214;
  
  fHist_Vx_ArrayFinder = new TH1F("fHist_Vx_ArrayFinder","",vxBin,VxCut[0],VxCut[1]);
  fHist_Vy_ArrayFinder = new TH1F("fHist_Vy_ArrayFinder","",vyBin,VyCut[0],VyCut[1]);
  fHist_Vz_ArrayFinder = new TH1F("fHist_Vz_ArrayFinder","",vzBin,   -10,    10);
  

  const Double_t dCent[20] = {0,0.5,1,2,3,4,5,6,7,8,10,15,20,30,40,50,60,70,80,90}; //not used in this code.

  for(int i=0;i<89;i++){
    fHist_znCx_V0_VxVy[i] = new TProfile2D(Form("fHist_znCx_V0_Run%d_Vz%d",runNums[i],1),"",NbinVt,0,NbinVt,vzBin,-10,10,"");
    fHist_znCy_V0_VxVy[i] = new TProfile2D(Form("fHist_znCy_V0_Run%d_Vz%d",runNums[i],1),"",NbinVt,0,NbinVt,vzBin,-10,10,"");
    fHist_znAx_V0_VxVy[i] = new TProfile2D(Form("fHist_znAx_V0_Run%d_Vz%d",runNums[i],1),"",NbinVt,0,NbinVt,vzBin,-10,10,"");
    fHist_znAy_V0_VxVy[i] = new TProfile2D(Form("fHist_znAy_V0_Run%d_Vz%d",runNums[i],1),"",NbinVt,0,NbinVt,vzBin,-10,10,"");
  }

for(int i=0;i<89;i++){ 
    fListHistos->Add(fHist_znCx_V0_VxVy[i]);
    fListHistos->Add(fHist_znCy_V0_VxVy[i]);
    fListHistos->Add(fHist_znAx_V0_VxVy[i]);
    fListHistos->Add(fHist_znAy_V0_VxVy[i]);
 }
 //-------------------------------------------------------------------
  
// Read the <Vx>, <Vy> and <V> for extreme cuts:
 fHist_Vx_vs_runnumber = (TProfile *) fZDCESEList->FindObject("fHist_Vx_vs_runnum");
 fHist_Vy_vs_runnumber = (TProfile *) fZDCESEList->FindObject("fHist_Vy_vs_runnum");

     if(!fHist_Vx_vs_runnumber || !fHist_Vy_vs_runnumber){
       printf("\n\n ****...<Vx>,<Vy> vs runnumber Histograms not found. Quit...***!!! \n\n");
       exit(0);
     }


  PostData(1,fListHistos); //posting in slot #1,
  printf("\n ... AliAnalysisTaskBernie::UserCreateOutPutObject() is called successfully... NbinVt = %d \n",NbinVt);

}

//________________________________________________________________________
void AliAnalysisTaskBernie::UserExec(Option_t *)
{
  //printf("\n ... AliAnalysisTaskBernie::UserExec() is being called 1...  \n");

  fHist_Event_count->Fill(1);
  
  /*if(!checkOnce){
  std::cout<<" runs = "<<std::endl;
  for(int i=0;i<89;i++)
  std::cout<<runNums[i]<<",";
  std::cout<<std::endl;
  checkOnce++;
  }*/
  
  AliAODEvent *aod = dynamic_cast<AliAODEvent*>(InputEvent());	   
  fEvent           = dynamic_cast<AliFlowEventSimple*>(GetInputData(1));

  fHist_Event_count->Fill(2);


  if(!aod)            return;
  fHist_Event_count->Fill(3);
 //printf("\n ... AliAnalysisTaskBernie:: no aod found.....  \n");


  //===========Rihan=========> get runindex:
  Int_t runNumber = aod->GetRunNumber();
  Int_t runindex = -1;

  for(int i=0;i<frunflag;i++){
   if(runNumber==runNums[i])
    {
      runindex = i;
      break;
    }
  }
  

 if(runNumber!=oldRunNum){
   maxVx_low  = fHist_Vx_vs_runnumber->GetBinContent(runindex+1) - 3.0*fHist_Vx_vs_runnumber->GetBinError(runindex+1)*sqrt(fHist_Vx_vs_runnumber->GetBinEntries(runindex+1));
   maxVx_high = fHist_Vx_vs_runnumber->GetBinContent(runindex+1) + 3.0*fHist_Vx_vs_runnumber->GetBinError(runindex+1)*sqrt(fHist_Vx_vs_runnumber->GetBinEntries(runindex+1));

   maxVy_low  = fHist_Vy_vs_runnumber->GetBinContent(runindex+1) - 3.0*fHist_Vy_vs_runnumber->GetBinError(runindex+1)*sqrt(fHist_Vy_vs_runnumber->GetBinEntries(runindex+1));
   maxVy_high = fHist_Vy_vs_runnumber->GetBinContent(runindex+1) + 3.0*fHist_Vy_vs_runnumber->GetBinError(runindex+1)*sqrt(fHist_Vy_vs_runnumber->GetBinEntries(runindex+1)); 

 oldRunNum = runNumber;
 printf("\n ... AliAnalysisTaskBernie::UserExec() New run number found = %d, index = %d  \n\n",runNumber,runindex);
}

  //fHist_EventperRun->Fill(runindex);
  
  AliAODVertex *pVertex = aod->GetPrimaryVertex();
  Double_t Vxyz[3]    =         {0,0,0};
  Vxyz[0]             = pVertex->GetX();
  Vxyz[1]             = pVertex->GetY();
  Vxyz[2]             = pVertex->GetZ();


  //------- Apply Necessary event cuts ---------
  if(Vxyz[2] >= 9.2    ||  Vxyz[2] <= -10.0)        return;   // ****** 9.2 cm Vz cut due to Correlation with centrality.
  fHist_Event_count->Fill(4);
  
  if(Vxyz[0] >= maxVx_high ||  Vxyz[0] <= maxVx_low) return;
  fHist_Event_count->Fill(5);
  
  if(Vxyz[1] >= maxVy_high  || Vxyz[1] <= maxVy_low) return;
  fHist_Event_count->Fill(6);
	    
 	   
  AliAODZDC *aodZDC = aod->GetZDCData();

  Float_t                            fZDCGainAlpha = 0.500;   // ********** fZDCgain alpha = 0.50 instead of 0.35 *********
  Float_t energyZNC  = (Float_t)  (aodZDC->GetZNCEnergy());
  Float_t energyZPC  = (Float_t)  (aodZDC->GetZPCEnergy());
  Float_t energyZNA  = (Float_t)  (aodZDC->GetZNAEnergy());
  Float_t energyZPA  = (Float_t)  (aodZDC->GetZPAEnergy());
  Float_t energyZEM1 = (Float_t) (aodZDC->GetZEM1Energy());
  Float_t energyZEM2 = (Float_t) (aodZDC->GetZEM2Energy());

  const Double_t * towZNC   = aodZDC->GetZNCTowerEnergy();
  const Double_t * towZPC   = aodZDC->GetZPCTowerEnergy();
  const Double_t * towZNA   = aodZDC->GetZNATowerEnergy();
  const Double_t * towZPA   = aodZDC->GetZPATowerEnergy();

  const Double_t * towZNClg = aodZDC->GetZNCTowerEnergyLR(); // Low gain something, should not be used.
  const Double_t * towZNAlg = aodZDC->GetZNATowerEnergyLR();

  Double_t towZPClg[5] = {0.,};
  Double_t towZPAlg[5] = {0.,};

  for(Int_t it=0; it<5; it++) {
    towZPClg[it] = 8*towZPC[it];
    towZPAlg[it] = 8*towZNA[it];
  }

  //***** Get centroid from ZDCs***************                                                               
  Double_t xyZNC[2]={999.,999.};
  Double_t xyZNA[2]={999.,999.};

  Float_t zncEnergy=0., znaEnergy=0.;

  for(Int_t i=0; i<5; i++){
    zncEnergy += towZNC[i];
    znaEnergy += towZNA[i];
  }

  //fEvent->SetZNCEnergy(towZNC[0]);
  //fEvent->SetZNAEnergy(towZNA[0]);

  Double_t AvTowerGain[8] = {1., 1., 1., 1., 1., 1., 1., 1.};

/*---------------------------------------------------------------    
  Int_t CenBin = GetCenBin(centrperc);
  Double_t zvtxpos[3]={0.,0.,0.};
  fFlowEvent->GetVertexPosition(zvtxpos);
  Int_t RunNum=fFlowEvent->GetRun();
  if(fTowerEqList) {
   if(RunNum!=fCachedRunNum) {
      for(Int_t i=0; i<8; i++) {
      fTowerGainEq[i] = (TH1D*)(fTowerEqList->FindObject(Form("Run %d",RunNum))->FindObject(Form("fhnTowerGainEqFactor[%d][%d]",RunNum,i)));
     }
   }
 }   
  Bool_t fUseMCCen = kFALSE; //rihan:hardcoded
  if (fUseMCCen) {
   if(aod->GetRunNumber() < 209122)
    aodZDC->GetZNCentroidInPbPb(1380., xyZNC, xyZNA);
    else
    aodZDC->GetZNCentroidInPbPb(2510., xyZNC, xyZNA);
  }
  else {
  //set tower gain equalization, if available     
  if(fTowerEqList) {
   for(Int_t i=0; i<8; i++)
   {
    if(fTowerGainEq[i])
    AvTowerGain[i] = fTowerGainEq[i]->GetBinContent(fTowerGainEq[i]->FindBin(centrperc));
   }
 }//---------------------------------------------------------------  */

  
  
  const Float_t x[4] = {-1.75,  1.75,-1.75, 1.75};
  const Float_t y[4] = {-1.75, -1.75, 1.75, 1.75};

  Float_t numXZNC=0., numYZNC=0., denZNC=0., wZNC;
  Float_t numXZNA=0., numYZNA=0., denZNA=0., wZNA;

  for(Int_t i=0; i<4; i++)
   {
    if(towZNC[i+1]>0.)
      {
       wZNC = TMath::Power(towZNC[i+1], fZDCGainAlpha)*AvTowerGain[i];
       numXZNC += x[i]*wZNC;
       numYZNC += y[i]*wZNC;
       denZNC  += wZNC;
       }

    if(towZNA[i+1]>0.) {
       wZNA = TMath::Power(towZNA[i+1], fZDCGainAlpha)*AvTowerGain[i+4];
       numXZNA += x[i]*wZNA;  
       numYZNA += y[i]*wZNA;
       denZNA  += wZNA;
       }
    }

  if(denZNC!=0) {
    xyZNC[0] = numXZNC/denZNC;
    xyZNC[1] = numYZNC/denZNC;
  }
   else{
     xyZNC[0]  = 999.;
     xyZNC[1]  = 999.;
     zncEnergy =   0.;
    }
  if(denZNA!=0) {
     xyZNA[0] = numXZNA/denZNA;
     xyZNA[1] = numYZNA/denZNA;
    }
   else{
     xyZNA[0]  = 999.;
     xyZNA[1]  = 999.;
     znaEnergy =   0.;
   }


  //----- Important: zdcA_X = -zdcA_X ---------
  xyZNA[0] = -1.*xyZNA[0];
  //===========================================
  

  if(sqrt(xyZNC[0]*xyZNC[0] + xyZNC[1]*xyZNC[1])>1.5)
   return;
  fHist_Event_count->Fill(7);

  if(sqrt(xyZNA[0]*xyZNA[0] + xyZNA[1]*xyZNA[1])>1.5) 
   return;
  fHist_Event_count->Fill(8);

  //this next cut should be removed
  //if(xyZNA[0]<1.e-5 && xyZNA[1]<1.e-5 && xyZNC[0]<1.e-5 && xyZNC[1]<1.e-5) return;

  fHist_Event_count->Fill(9);

  
  fHist_ZDCn_A_XY->Fill(xyZNA[0],xyZNA[1]);
  fHist_ZDCn_C_XY->Fill(xyZNC[0],xyZNC[1]);
 

  
/*// ----------------  ZDC multiplicities -----------------------
//              posted at the ende of this code
//-------------------- ZDC Mutliplicities -----------------------   */



 //----- calculate RefMult for event ---------
    AliFlowTrackSimple*   pTrack = NULL;
    Int_t iTracks = fEvent->NumberOfTracks();

    Double_t    Qnx_TPC[3]  = {0,};
    Double_t    Qny_TPC[3]  = {0,};
    Double_t         dPhi,dPt,dEta;

    Int_t             nRefMult = 0;
    Int_t             npoiMult = 0;

    
    for(Int_t i=0; i<iTracks; i++)
    {
    pTrack     =  fEvent->GetTrack(i);
    if (!pTrack)             continue;
    
    dPhi      =         pTrack->Phi();
    dPt       =         pTrack-> Pt();
    dEta      =         pTrack->Eta();
  //dCharge   =      pTrack->Charge();
    
    if(dPt<0.2 || dPt>10.0)  continue;
    if(fabs(dEta)>0.8)       continue;

    fHistQA_etaphi->Fill(dPhi,dEta);

    Qnx_TPC[0] += TMath::Cos(2.*dPhi);
    Qny_TPC[0] += TMath::Sin(2.*dPhi); 
        
    npoiMult++;
    

    if(fabs(dEta)<0.5)
    nRefMult++;
    }//track loop ends


    double psi2 = 0.5*TMath::ATan2(Qny_TPC[0],Qnx_TPC[0]);
    if(psi2<0) psi2 += TMath::Pi();

    fHist_EventPlane2->Fill(psi2);

     psi2 = TMath::ATan2(xyZNC[1],xyZNC[0]);
    if(psi2<0) psi2 += 2.*TMath::Pi();
    fHist_Psi1_zdnC->Fill(psi2);

     psi2 = TMath::ATan2(xyZNA[1],xyZNA[0]);
    if(psi2<0) psi2 += 2.*TMath::Pi();
    fHist_Psi1_zdnA->Fill(psi2);

    
   Double_t      EvtCent = fEvent->GetCentrality();
   fHist_Centr_count     ->          Fill(EvtCent);
    
 //fill q vectors for TPC recentering
   //fHist_QnxRecent->Fill(EvtCent,runindex,(Qnx_TPC[0]/(double)npoiMult));
   //fHist_QnyRecent->Fill(EvtCent,runindex,(Qny_TPC[0]/(double)npoiMult));
       
    //---------------------------------


  Int_t nTracks = aod->GetNumberOfTracks();

  fHist_Vertex_Vz->Fill(Vxyz[2]);
  fHist_Vertex_XY->Fill(Vxyz[0],Vxyz[1]);

  
  fHist_Vx_vs_runnum->Fill(runindex,Vxyz[0]);
  fHist_Vy_vs_runnum->Fill(runindex,Vxyz[1]);
  fHist_Vz_vs_runnum->Fill(runindex,Vxyz[2]);
  
  fHist_tracks_vs_runnum->Fill(runindex,nTracks);

    Int_t indexVx = fHist_Vx_ArrayFinder->FindBin(Vxyz[0]);
    Int_t indexVy = fHist_Vy_ArrayFinder->FindBin(Vxyz[1]);
    Int_t indexVz = fHist_Vz_ArrayFinder->FindBin(Vxyz[2]);

    //int tVertexBin = (indexVz-1)*20*20 + (indexVy-1)*20 + indexVx - 1; //
    Double_t tVertexBin =  (Double_t) (indexVy-1)*20 + indexVx - 0.5     ; // use double because fillin in a bin

    fHist_vBincount->Fill(tVertexBin);


    indexVz = indexVz - 1;

    //Int_t iCent = (int) EvtCent;

  fHist_znCx_V0_VxVy[runindex]->Fill(tVertexBin,Vxyz[2],xyZNC[0]);
  fHist_znCy_V0_VxVy[runindex]->Fill(tVertexBin,Vxyz[2],xyZNC[1]);
  fHist_znAx_V0_VxVy[runindex]->Fill(tVertexBin,Vxyz[2],xyZNA[0]);
  fHist_znAy_V0_VxVy[runindex]->Fill(tVertexBin,Vxyz[2],xyZNA[1]);
  

  //printf("\n ... AliAnalysisTaskBernie::UserExec() before printing Vycut = %f  \n",maxVy_high);

  if(fievent%50==0)
    std::cout<<fievent<<" cent = "<<EvtCent<<" iRun = "<<runindex<<" lowVx = "<<maxVx_low<<"\tVx = "<<Vxyz[0]<<"\thigVx = "
             <<maxVx_high<<"\tlowVy = "<<maxVy_low<<"\t Vy = "<<Vxyz[1]<<"\t higVy = "<<maxVy_high<<std::endl;

  PostData(1,fListHistos);
    
  fievent++;

  //printf("\n ... AliAnalysisTaskBernie::UserExec() is being called 5...  \n");
} 

//________________________________________________________________________
void AliAnalysisTaskBernie::Terminate(Option_t *)
{
  // Called once at the end of the query
       /* AliFlowAnalysisIDCSP *fSPTerm = new AliFlowAnalysisIDCSP();
	fListHistos = (TList*) GetOutputData(1);
        if(fListHistos)
	  {
	   fSPTerm->GetOutputHistograms(fListHistos);
	   fSPTerm->Finish();
	   PostData(1,fListHistos);
	  } 
        else
	{ 
          std::cout << "histgram list pointer is empty in Scalar Product" << endl;
         } */

  //delete fListHistos;
	printf("\n ... AliAnalysisTaskBernie::Terminate() is being called ...  \n");
}

/*// ----------------  ZDC multiplicities -----------------------
  Float_t MulA=0., MulC=0.;
  for(Int_t i=0; i<4; i++) {
    if(towZNC[i+1]>0.) {
	MulC += TMath::Power(towZNC[i+1], fZDCGainAlpha)*AvTowerGain[i+4];
    }
    if(towZNA[i+1]>0.) {
	MulA += TMath::Power(towZNA[i+1], fZDCGainAlpha)*AvTowerGain[i+4];
    }
  }
//fhZNCcentroid->Fill(xyZNC[0], xyZNC[1]);
//fhZNAcentroid->Fill(xyZNA[0], xyZNA[1]);
  Double_t testC[2] = {1,0};
  Double_t testA[2] = {1,0};
//fEvent->SetZDC2Qsub(xyZNC,MulC,xyZNA,MulA);
  Float_t tdcSum  = aodZDC->GetZDCTimeSum();
  Float_t tdcDiff = aodZDC->GetZDCTimeDiff();
  Double_t asymmetry = -999.;
  if((energyZNC+energyZNA)>0.)
  asymmetry = (energyZNC-energyZNA)/(energyZNC+energyZNA);
//----------------------- ZDC Mutliplicities -----------------------   */
