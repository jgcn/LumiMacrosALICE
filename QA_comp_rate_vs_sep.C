#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// Macro to produce quality-assurance (QA) plots related to
// the rate
//-------------------------------------------------------

void QA_comp_rate_vs_sep(Int_t Fill, const char *rate_name1, const char *rate_type1,
		            const char *sep_type1, const char *rate_name2, const char *rate_type2,
		            const char *sep_type2,
			 Int_t scan_type, Int_t scan, Int_t bc)
// scan_type: 1 => x-scan; 2 => y-scan
{
	// initialize
	Set_input_file_names(Fill);
	Set_pointers_to_input_files_and_trees();

	//-------------------------------------------------------
	// first rate
	//-------------------------------------------------------
	
	// --> rate file
	char *file_name = new char[kg_string_size];
	//	TString file_name;
	if (scan_type == 1) sprintf(file_name,"../Fill-%d/%sRate_%s_x_Scan_%d.root",g_vdm_Fill,rate_type1,rate_name1,scan);
        if (scan_type == 2) sprintf(file_name,"../Fill-%d/%sRate_%s_y_Scan_%d.root",g_vdm_Fill,rate_type1,rate_name1,scan);
	TFile *rate_file1 = new TFile(file_name);
	// --> separation file
	if (scan_type == 1) sprintf(file_name,"../Fill-%d/%sSep_x_Scan_%d.root",g_vdm_Fill,sep_type1,scan);
	if (scan_type == 2) sprintf(file_name,"../Fill-%d/%sSep_y_Scan_%d.root",g_vdm_Fill,sep_type1,scan);
	TFile *sep_file1 = new TFile(file_name);
	// --> get the trees
	TTree *rate_tree1 = (TTree *) rate_file1->Get("Rate");
	TTree *sep_tree1 = (TTree *) sep_file1->Get("Separations");

	// Next step: prepare variables to be fitted
	// --> get number of separations
	Int_t n_sep = FindNumberSeparations(scan_type, scan);
	// --> reserve space for rates and separations
	Double_t *rate1 = new Double_t[n_sep];
	Double_t *rate_error1 = new Double_t[n_sep];  
	Double_t *sep1 = new Double_t[n_sep];
	// --> set branches
	rate_tree1->ResetBranchAddresses();
	rate_tree1->SetBranchAddress("rate",rate1);
	rate_tree1->SetBranchAddress("rate_error",rate_error1);
	sep_tree1->ResetBranchAddresses();
	sep_tree1->SetBranchAddress("separation",sep1);

	// get the right bunch crossing
	rate_tree1->GetEntry(bc);
	sep_tree1->GetEntry(bc);

	// fill graph
	TGraphErrors *gr1 = new TGraphErrors(n_sep,sep1,rate1,NULL,rate_error1);
	gr1->SetMarkerStyle(20);
	gr1->SetMarkerColor(kBlue);
	//-------------------------------------------------------
	// second rate
	//-------------------------------------------------------
	
	// --> rate file
	if (scan_type == 1) sprintf(file_name,"../Fill-%d/%sRate_%s_x_Scan_%d.root",g_vdm_Fill,rate_type2,rate_name2,scan);
	if (scan_type == 2) sprintf(file_name,"../Fill-%d/%sRate_%s_y_Scan_%d.root",g_vdm_Fill,rate_type2,rate_name2,scan);
	TFile *rate_file2 = new TFile(file_name);
	// --> separation file
	if (scan_type == 1) sprintf(file_name,"../Fill-%d/%sSep_x_Scan_%d.root",g_vdm_Fill,sep_type2,scan);
	if (scan_type == 2) sprintf(file_name,"../Fill-%d/%sSep_y_Scan_%d.root",g_vdm_Fill,sep_type2,scan);
	TFile *sep_file2 = new TFile(file_name);
	// --> get the trees
	TTree *rate_tree2 = (TTree *) rate_file2->Get("Rate");
	TTree *sep_tree2 = (TTree *) sep_file2->Get("Separations");

	// Next step: prepare variables to be fitted
	// --> reserve space for rates and separations
	Double_t *rate2 = new Double_t[n_sep];
	Double_t *rate_error2 = new Double_t[n_sep];  
	Double_t *sep2 = new Double_t[n_sep];
	// --> set branches
	rate_tree2->ResetBranchAddresses();
	rate_tree2->SetBranchAddress("rate",rate2);
	rate_tree2->SetBranchAddress("rate_error",rate_error2);
	sep_tree2->ResetBranchAddresses();
	sep_tree2->SetBranchAddress("separation",sep2);

	// get the right bunch crossing
	rate_tree2->GetEntry(bc);
	sep_tree2->GetEntry(bc);

	// fill graph
	TGraphErrors *gr2 = new TGraphErrors(n_sep,sep2,rate2,NULL,rate_error2);
	gr2->SetMarkerStyle(20);
	gr2->SetMarkerColor(kRed);

	//-------------------------------------------------------
	// ratio
	//-------------------------------------------------------

	Double_t *ratio = new Double_t[n_sep];
	Double_t *ratio_err = new Double_t[n_sep];
	int mid = (n_sep>>1)+1;
	double r01 = rate1[mid];
	double r02 = rate2[mid];	
	for(int i=0;i<n_sep;i++) {
	  ratio[i] = (rate1[i]/r01)/(rate2[i]/r02);
	  // error computation not correct
	  ratio_err[i] =ratio[i]*TMath::Sqrt(TMath::Power(rate_error1[i]/rate1[i],2)+TMath::Power(rate_error2[i]/rate2[i],2));
	}
	TGraphErrors *grR = new TGraphErrors(n_sep,sep2,ratio,NULL,ratio_err);
	grR->SetMarkerStyle(20);
	grR->SetMarkerColor(kBlack);
	grR->Draw();
	//-------------------------------------------------------
	// now plot
	//-------------------------------------------------------
	

	// define the limits for the plot
	// --> separation
	Double_t sep_min = 0;
	Double_t sep_max = 0;
	for(Int_t i=0;i<n_sep;i++) {
		if(sep1[i]<sep_min) sep_min=sep1[i];
		if(sep1[i]>sep_max) sep_max=sep1[i];
		if(sep2[i]<sep_min) sep_min=sep2[i];
		if(sep2[i]>sep_max) sep_max=sep2[i];
	}
	sep_min*=1.2;
	sep_max*=1.2;  
	// --> find maximum rate
	Double_t rate_max = 0;
	for(Int_t i=0;i<n_sep;i++) {
		if ((rate1[i]+rate_error1[i])>rate_max) rate_max = rate1[i]+rate_error1[i];
		if ((rate2[i]+rate_error2[i])>rate_max) rate_max = rate2[i]+rate_error2[i];

	}
	rate_max *= 1.2; //consider 20% larger limit

	// plot TGraph
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	TCanvas *rvss_C = new TCanvas("rvss_C","rate versus separation",900,400);
	rvss_C->Divide(2,1);
	rvss_C->cd(1);
	TH1F* frame = gPad->DrawFrame(sep_min,0.001,sep_max,rate_max);
	frame->SetTitle(";separation (mm); rate (Hz)");
	gr1->Draw("p,e1,same");
	gr2->Draw("p,e1,same");
	TLatex* txt = new TLatex();
	txt->SetTextAlign(11);
	txt->SetNDC();
	txt->SetTextSize(0.035);
	txt->SetTextFont(42);
	txt->SetTextColor(kBlue);
	txt->DrawLatex(0.25,0.8,rate_name1);
	txt->SetTextColor(kRed);
	txt->DrawLatex(0.25,0.75,rate_name2);
   
	rvss_C->cd(2);
	grR->Draw("");
	TH1 *h = (TH1*) grR->GetHistogram();
	h->SetTitle(Form(";separation (mm);((%s/%s(0))/(%s/%s(0))",rate_name1,rate_name1,rate_name2,rate_name2));
	//	rvss_C->Print(Form("c2a_QA_comp_rate_Fill%i_%s_%s_%s_scanT%i_scan%i_bc%i.%s",
	//			Fill, rate_name, rate_type, sep_type, scan_type, scan, bc, FFormat));
	
}
