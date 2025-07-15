#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// Macro to produce quality-assurance (QA) plots related to
// the rate
//-------------------------------------------------------

void test(int fill, const char *rate_name) {
  // get tree
  TFile *f = new TFile(Form("../Fill-%d/vdm_time_9128_DUMMY_1_v3.root",fill));
  TTree *t = (TTree *) f->Get("VdM");
  t->ResetBranchAddresses();
  int orbit;
  int aqflag;
  double timerel;
  double *counter  = new Double_t[3564];
  t->SetBranchAddress("orbits",&orbit);
  t->SetBranchAddress("aqflag",&aqflag);
  t->SetBranchAddress("timerel",&timerel);
  t->SetBranchAddress(Form("%scounts",rate_name),counter);

  // define histogram
  int n = t->GetEntries();
  TH2D *h = new TH2D("countsTimerel",";timerel;bc;counts",n,-0.5,n-0.5,3564,-0.5,3564-0.5);
  // loop over tree
  for(int i=0;i<n;i++) {
    t->GetEntry(i);
    if (aqflag==0) continue;
    for(int j=0;j<3564;j++) h->SetBinContent(i+1,j+1,counter[j]);
  }

  // plot
  h->Draw("zcol");
}


void QA_get_ACE_rate(int fill, const char *rate_name, int scan_type, int scan)
// right now several things are fixed
{
  // initialize
  Set_input_file_names(fill);
  Set_pointers_to_input_files_and_trees();

  // --> rate file
  char *file_name = new char[kg_string_size];
  //	TString file_name;
  if (scan_type == 1) sprintf(file_name,"../Fill-%d/RawRate_%s_x_Scan_%d_allBC.root",fill,rate_name,scan);
  if (scan_type == 2) sprintf(file_name,"../Fill-%d/RawRate_%s_y_Scan_%d_allBC.root",fill,rate_name,scan);
  std::cout << file_name << " " << fill << " " << scan << " " << rate_name << std::endl;
  TFile *rate_file = new TFile(file_name);
  // --> separation file
  if (scan_type == 1) sprintf(file_name,"../Fill-%d/NomSep_x_Scan_%d.root",fill,scan);
  if (scan_type == 2) sprintf(file_name,"../Fill-%d/NomSep_y_Scan_%d.root",fill,scan);
  std::cout << file_name << " " << fill << " " << scan << " " << rate_name << std::endl;
  TFile *sep_file = new TFile(file_name);
  
  // --> get the trees
  TTree *rate_tree = (TTree *) rate_file->Get("Rate");
  TTree *sep_tree = (TTree *) sep_file->Get("Separations");
  std::cout << " rate tree has " << rate_tree->GetEntries()  << " entries" << std::endl;
  std::cout << " sep tree has " << sep_tree->GetEntries()  << " entries" << std::endl;
  
  // Next step: prepare variables to be fitted
  // --> get number of separations
  Int_t n_sep = FindNumberSeparations(scan_type, scan);
  std::cout << " There are " << n_sep << " separtions " << std::endl;
  
  // --> reserve space for rates and separations
  Double_t *rate = new Double_t[n_sep];
  Double_t *sep = new Double_t[n_sep];
  // --> set branches
  rate_tree->ResetBranchAddresses();
  rate_tree->SetBranchAddress("rate",rate);
  sep_tree->ResetBranchAddresses();
  sep_tree->SetBranchAddress("separation",sep);
  sep_tree->GetEntry(0);
  
  // get the fill scheme 
  int fillScheme[3564];
  ifstream ifs;
  ifs.open(Form("fillScheme%d.txt",fill));
  for(int i=0;i<3564;i++) ifs >> fillScheme[i];
  ifs.close();
  
  // set up storage space
  Double_t *rateA = new Double_t[n_sep];
  Double_t *rateB = new Double_t[n_sep];  
  Double_t *rateC = new Double_t[n_sep];	
  Double_t *rateE = new Double_t[n_sep];
  int nA=0;
  int nB=0;
  int nC=0;
  int nE=0;
  for(int j=0;j<n_sep;j++) rateA[j]=rateB[j]=rateC[j]=rateE[j]=0;
  
  // get the rates
  for(int i=0;i<3564;i++) {
    rate_tree->GetEntry(i);
    
    for(int j=0;j<n_sep;j++) {
      if (fillScheme[i]==0) rateA[j]+=rate[j];
      if (fillScheme[i]==1) rateB[j]+=rate[j];      
      if (fillScheme[i]==2) rateC[j]+=rate[j];
      if (fillScheme[i]==3) rateE[j]+=rate[j];	    
    }
    if (fillScheme[i]==0) nA++;
    if (fillScheme[i]==1) nB++;    
    if (fillScheme[i]==2) nC++;
    if (fillScheme[i]==3) nE++;	    
  }
  // normalise the rates
  for(int j=0;j<n_sep;j++) {
    std::cout << j << " " << rateA[j] << " " << rateB[j] << " "<< rateC[j] << " " << rateE[j] << std::endl; 
    rateA[j]/=((float) nA);
    rateB[j]/=((float) nB);
    rateC[j]/=((float) nC);
    rateE[j]/=((float) nE);
    std::cout << " == " << j << " " << rateA[j] << " " << rateB[j] << " "<< rateC[j] << " " << rateE[j] << std::endl; 
  }
  std::cout  << " nA " << nA << " nB " << nB  << " nC " << nC << " nE " << nE << std::endl; 

  
  // fill graphs
  TGraph *grA = new TGraph(n_sep,sep,rateA);
  grA->SetMarkerStyle(20);
  grA->SetMarkerColor(kBlue);
  TGraph *grB = new TGraph(n_sep,sep,rateB);
  grB->SetMarkerStyle(24);
  grB->SetMarkerSize(1.5);  
  grB->SetMarkerColor(kMagenta);
  TGraph *grC = new TGraph(n_sep,sep,rateC);
  grC->SetMarkerStyle(20);
  grC->SetMarkerColor(kRed);
  TGraph *grE = new TGraph(n_sep,sep,rateE);
  grE->SetMarkerStyle(20);
  grE->SetMarkerColor(kBlack);
  
  //-------------------------------------------------------
  // now plot
  //-------------------------------------------------------
  
  // define the limits for the plot
  // --> separation
  Double_t sep_min = 0;
  Double_t sep_max = 0;
  for(Int_t i=0;i<n_sep;i++) {
    if(sep[i]<sep_min) sep_min=sep[i];
    if(sep[i]>sep_max) sep_max=sep[i];
  }
  sep_min*=1.2;
  sep_max*=1.2;  
  // --> find maximum rate
  Double_t rate_max = 0;
  for(Int_t i=0;i<n_sep;i++) {
    if (rateA[i]>rate_max) rate_max = rateA[i];
    if (rateB[i]>rate_max) rate_max = rateB[i];
    if (rateC[i]>rate_max) rate_max = rateC[i];
    if (rateE[i]>rate_max) rate_max = rateE[i];		    
  }
  rate_max *= 1.2; //consider 20% larger limit

  // plot TGraph
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas *rvss_C = new TCanvas(Form("rates%s",rate_name),Form("%s rate ACE",rate_name),900,600);
  TH1F* frame = gPad->DrawFrame(sep_min,0.001,sep_max,rate_max);
  frame->SetTitle(";separation (mm); rate (Hz)");
  grA->Draw("p,same");
  grB->Draw("p,same");  
  grC->Draw("p,same");
  grE->Draw("p,same");
  
  // clean up
  delete [] rate;
  delete [] sep;
  delete [] rateA; 
  delete [] rateC;
  delete [] rateE;  
}
