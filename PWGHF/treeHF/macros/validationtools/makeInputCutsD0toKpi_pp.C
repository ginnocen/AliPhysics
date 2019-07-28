#include <Riostream.h>
#include <TFile.h>
#include <AliRDHFCutsD0toKpi.h>
#include <TClonesArray.h>
#include <TParameter.h>


//setting my cut values
//cuts order
//     printf("    |M-MD0| [GeV]    < %f\n",fD0toKpiCuts[0]);
//     printf("    dca    [cm]  < %f\n",fD0toKpiCuts[1]);
//     printf("    cosThetaStar     < %f\n",fD0toKpiCuts[2]);
//     printf("    pTK     [GeV/c]    > %f\n",fD0toKpiCuts[3]);
//     printf("    pTpi    [GeV/c]    > %f\n",fD0toKpiCuts[4]);
//     printf("    |d0K|  [cm]  < %f\n",fD0toKpiCuts[5]);
//     printf("    |d0pi| [cm]  < %f\n",fD0toKpiCuts[6]);
//     printf("    d0d0  [cm^2] < %f\n",fD0toKpiCuts[7]);
//     printf("    cosThetaPoint    > %f\n",fD0toKpiCuts[8]);
//     printf("    |cosThetaPointXY| < %f\n",fD0toKpiCuts[9]);
//     printf("    NormDecayLenghtXY    > %f\n",fD0toKpiCuts[10]);

AliRDHFCutsD0toKpi* makeInputCutsD0toKpi_pp(TString nameCuts="D0toKpiValCuts")
{
  AliRDHFCutsD0toKpi* cutsD0toKpi=new AliRDHFCutsD0toKpi();
  cutsD0toKpi->SetName(nameCuts.Data());
  cutsD0toKpi->SetTitle(nameCuts.Data());
  
  //single track cuts
  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetPtRange(0.3,1.e10);
  cutsD0toKpi->AddTrackCuts(esdTrackCuts);
  cutsD0toKpi->SetSelectCandTrackSPDFirst(kTRUE, 3);
  
  const Int_t nvars=11;
  
  const Int_t nptbins=1;
  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  ptbins[0]=1.;
  ptbins[1]=10.;
  
  //m    dca      cost*  ptk ptpi  d0k          d0pi       d0d0          cosp  cosxy normdxy
  Float_t cutsMatrixD0toKpiStand[nptbins][nvars]=  
    {{0.400,300.*1E-4, 0.8, 0.7, 0.7, 1000.*1E-4, 1000.*1E-4, -8000. *1E-8, 0.85,  0.,0.}};
  
  //CREATE TRANSPOSE MATRIX...REVERSE INDICES as required by AliRDHFCuts
  Float_t **cutsMatrixTransposeStand=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++)cutsMatrixTransposeStand[iv]=new Float_t[nptbins];
  
  for (Int_t ibin=0;ibin<nptbins;ibin++){
    for (Int_t ivar = 0; ivar<nvars; ivar++){
      cutsMatrixTransposeStand[ivar][ibin]=cutsMatrixD0toKpiStand[ibin][ivar];
    }
  }
  
  cutsD0toKpi->SetPtBins(nptbins+1,ptbins);
  cutsD0toKpi->SetGlobalIndex(nvars,nptbins);
  cutsD0toKpi->SetMinPtCandidate(1.);
  cutsD0toKpi->SetCuts(nvars,nptbins,cutsMatrixTransposeStand);
  
  cutsD0toKpi->SetUseSpecialCuts(kTRUE);
  
  for(Int_t iv=0;iv<nvars;iv++) delete [] cutsMatrixTransposeStand[iv];
  delete [] cutsMatrixTransposeStand;
  cutsMatrixTransposeStand=NULL;
  
  //pid settings
  AliAODPidHF* pidObj=new AliAODPidHF();
  //pidObj->SetName("pid4D0");
  Int_t mode=1;
  const Int_t nlims=2;
  Double_t plims[nlims]={0.6,0.8}; //TPC limits in momentum [GeV/c]
  Bool_t compat=kTRUE; //effective only for this mode
  Bool_t asym=kTRUE;
  Double_t sigmas[5]={2.,1.,0.,3.,0.}; //to be checked and to be modified with new implementation of setters by Rossella
  pidObj->SetAsym(asym);// if you want to use the asymmetric bands in TPC
  pidObj->SetMatch(mode);
  pidObj->SetPLimit(plims,nlims);
  pidObj->SetSigma(sigmas);
  pidObj->SetCompat(compat);
  pidObj->SetTPC(kTRUE);
  pidObj->SetTOF(kTRUE);
  pidObj->SetPCompatTOF(1.5);
  pidObj->SetSigmaForTPCCompat(3.);
  pidObj->SetSigmaForTOFCompat(3.);
  pidObj->SetOldPid(kFALSE);
  cutsD0toKpi->SetPidHF(pidObj);
  
  Bool_t pidflag=kTRUE;
  cutsD0toKpi->SetUsePID(pidflag);
  if(pidflag) cout<<"PID is used for analysis cuts"<<endl;
  else cout<<"PID is not used for analysis cuts"<<endl;
  
  cutsD0toKpi->SetRemoveDaughtersFromPrim(kTRUE);
  
  //event selection
  cutsD0toKpi->SetUsePhysicsSelection(kTRUE);
  cutsD0toKpi->SetTriggerClass("");
  cutsD0toKpi->SetTriggerMask(AliVEvent::kINT7);
  cutsD0toKpi->SetOptPileup(AliRDHFCuts::kRejectMVPileupEvent);
  cutsD0toKpi->SetMinContribPileupMV(5);
  cutsD0toKpi->SetMaxVtxZ(10.);
  cutsD0toKpi->SetCutOnzVertexSPD(3);
  cutsD0toKpi->SetMinVtxContr(1);
  
  cout<<"This is the object I'm going to save:"<<endl;
  cutsD0toKpi->SetName(nameCuts.Data());
  cutsD0toKpi->SetTitle(nameCuts.Data());
  cutsD0toKpi->PrintAll();
  return cutsD0toKpi;
}
void makefiles(){
  TFile* fout=new TFile("D0toKpiCutsValidation.root", "recreate");   //set this!!
  fout->cd();
  AliRDHFCutsD0toKpi* D0toKpiAnalysisCuts =(AliRDHFCutsD0toKpi*) makeInputCutsD0toKpi_pp("D0toKpiAnalysisCuts");
  AliRDHFCutsD0toKpi* D0toKpiFilteringCuts =(AliRDHFCutsD0toKpi*) makeInputCutsD0toKpi_pp("D0toKpiFilteringCuts");
  D0toKpiAnalysisCuts->Write();
  D0toKpiFilteringCuts->Write();
  fout->Close();
}
