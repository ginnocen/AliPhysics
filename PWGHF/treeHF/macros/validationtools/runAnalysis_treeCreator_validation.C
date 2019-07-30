#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"


#include "AliPhysicsSelectionTask.h"
#include "AliMultSelectionTask.h"
#include "AliAnalysisTaskPIDResponse.h"

#include "AliAnalysisTaskSEHFTreeCreator.h"
#endif

using namespace std;
enum system {kpp, kpPb, kPbPb};

//SETTINGS
//************************************
Bool_t runLocal=kTRUE;             
TString pathToLocalAODfiles;
TString runMode="full";          

TString aliPhysVersion="vAN-20190701_ROOT6-1";
TString gridDataDir;
TString gridDataPattern="/pass1/AOD194";

// Alien output directory
TString gridWorkingDir="testNtupleCreator";
TString gridOutputDir="output";

//run numbers
const Int_t nruns = 1;
Int_t runlist[nruns] = {246994};

//Task configuration
TString cutFile="D0toKpiCutsValidation.root";  

void runAnalysis_treeCreator_validation(Bool_t isRunOnMC=kFALSE)
{
  TGrid::Connect("alien://");
  
  Int_t System=-1;
  if (isRunOnMC == kTRUE) {
    pathToLocalAODfiles = "/data2/highmultppAOD/mc/LHC18f4a";
    gridDataDir="/alice/sim/2018/LHC18f4a";
    System=kpp;
    std::cout<<"Running on MC"<<std::endl;
  }
  if (isRunOnMC == kFALSE){
    pathToLocalAODfiles = "/data2/validationAOD/001";
    gridDataDir="/alice/data/2018/LHC18f";
    System=kpp;
    std::cout<<"Running on data"<<std::endl;
  }
  Bool_t local = runLocal;
  gInterpreter->ProcessLine(".include $ROOTSYS/include");
  gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
  gInterpreter->ProcessLine(".include $ALICE_PHYSICS/include");
  // create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskExample");
  AliAODInputHandler *aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);
  
  AliPhysicsSelectionTask *physseltask = reinterpret_cast<AliPhysicsSelectionTask *>(gInterpreter->ProcessLine(Form(".x %s (%d)", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"),isRunOnMC)));
  AliAnalysisTaskPIDResponse *pidResp = reinterpret_cast<AliAnalysisTaskPIDResponse *>(gInterpreter->ProcessLine(Form(".x %s (%d)", gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"),isRunOnMC)));
  
  if(System==kpPb || System==kPbPb) {
    AliMultSelectionTask *multSel = reinterpret_cast<AliMultSelectionTask *>(gInterpreter->ProcessLine(Form(".x %s", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"))));
    //multSel->SetAlternateOADBforEstimators("LHC15o-DefaultMC-HIJING");
    multSel->SetAlternateOADBFullManualBypassMC("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/data/OADB-LHC18q-DefaultMC-HIJING.root");
  }
  //AliMultSelectionTask *multSel = reinterpret_cast<AliMultSelectionTask *>(gInterpreter->ProcessLine(Form(".x %s", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C")))); 
  AliAnalysisTaskSED0Mass *taskD0 = reinterpret_cast<AliAnalysisTaskSED0Mass *>(gInterpreter->ProcessLine(Form(".x %s (%d,%d,%d,%d,%d,%d,%d,%d,\"%s\",\"%s\",\"%s\")",gSystem->ExpandPathName("$ALICE_PHYSICS/PWGHF/vertexingHF/macros/AddTaskD0Mass.C"),0,isRunOnMC,kFALSE,kFALSE,0,0,0,0,"_010",cutFile.Data(),"D0toKpiAnalysisCuts")));
 // AliAnalysisTaskSECleanupVertexingHF *taskclean = reinterpret_cast<AliAnalysisTaskSECleanupVertexingHF*>(gInterpreter->ProcessLine(Form(".x %s",gSystem->ExpandPathName("$ALICE_PHYSICS/PWGHF/vertexingHF/macros/AddTaskCleanupVertexingHF.C"))));

  AliAnalysisTaskSEHFTreeCreator *task = reinterpret_cast<AliAnalysisTaskSEHFTreeCreator*>(gInterpreter->ProcessLine(Form(".x %s (%d,%d,\"%s\",\"%s\", %d,%d,%d,%d,%d,%d,%d,%d,%d)",gSystem->ExpandPathName("$ALICE_PHYSICS/PWGHF/treeHF/macros/AddTaskHFTreeCreator.C"),isRunOnMC, 1, "HFTreeCreator", cutFile.Data(), 1, kTRUE, kTRUE,kTRUE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE)));
  if(System==kPbPb) {
    AliAnalysisTaskSECleanupVertexingHF *taskclean =reinterpret_cast<AliAnalysisTaskSECleanupVertexingHF *>(gInterpreter->ProcessLine(Form(".x %s", gSystem->ExpandPathName("$ALICE_PHYSICS/PWGHF/vertexingHF/macros/AddTaskCleanupVertexingHF.C"))));
  }
   
  if(!mgr->InitAnalysis()) return;
  mgr->SetDebugLevel(2);
  mgr->PrintStatus();
     
  // if you want to run locally, we need to define some input
  TChain* chainAOD = new TChain("aodTree");
  TChain *chainAODfriend = new TChain("aodTree");
  chainAOD->Add(Form("%s/AliAOD.root",pathToLocalAODfiles.Data()));
  chainAODfriend->Add(Form("%s/AliAOD.VertexingHF.root",pathToLocalAODfiles.Data()));
  chainAOD->AddFriend(chainAODfriend);
  mgr->SetUseProgressBar(1, 100000000);
  mgr->StartAnalysis("local", chainAOD);
}
