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
	// std::cout << " r02/r01 = " << (r02/r01) << endl;
	for(int i=0;i<n_sep;i++) {
	  // ratio[i] = (rate1[i]/r01)/(rate2[i]/r02);
	  ratio[i] = (rate2[i]>0?rate1[i]/rate2[i]:0);
	  // error computation not correct
	  ratio_err[i] =  (rate2[i]>0?
			   ratio[i]*TMath::Sqrt(TMath::Power(rate_error1[i]/rate1[i],2)+TMath::Power(rate_error2[i]/rate2[i],2))
			   :0);
	}
	TGraphErrors *grR = new TGraphErrors(n_sep,sep2,ratio,NULL,ratio_err);
	grR->SetMarkerStyle(20);
	grR->SetMarkerColor(kBlack);
	//	grR->Draw();
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
	gStyle->SetOptTitle(1);
	gStyle->SetOptFit(1);

	TCanvas *rvss_C = new TCanvas(Form("rate1_%s_rate2_%s_fill_%i_scan_%i_dir_%i_bc_%i",
					   rate_name1,rate_name2,Fill,scan,scan_type,bc),"rate versus separation",900,400);
	rvss_C->Divide(2,1);
	rvss_C->cd(1);
	TH1F* frame = gPad->DrawFrame(sep_min,0.001,sep_max,rate_max);
	frame->SetTitle(Form("fill %i, scan %i, dir %i, bc %i;separation (mm); rate (Hz)",Fill,scan,scan_type,bc));
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
	grR->Fit("pol0");
	TH1 *h = (TH1*) grR->GetHistogram();
	h->SetTitle(Form(";separation (mm);%s/%s",rate_name1,rate_name2));
       	rvss_C->Print(Form("aodPlots/rate1_%s_rate2_%s_fill_%i_bc_%i_scan_%i_dir_%i.pdf",
			   rate_name1,rate_name2,Fill,bc,scan,scan_type));
	
}


void doFill(const char *rate_name1, const char *rate_name2, int fill, int nBC = 20)
{
  for(int iBC=0;iBC<nBC;iBC++) {  
    QA_comp_rate_vs_sep(fill,rate_name1,"Raw","Nom",rate_name2,"Raw","Nom",1,0,iBC);
    QA_comp_rate_vs_sep(fill,rate_name1,"Raw","Nom",rate_name2,"Raw","Nom",2,0,iBC);
    QA_comp_rate_vs_sep(fill,rate_name1,"Raw","Nom",rate_name2,"Raw","Nom",1,1,iBC);
    QA_comp_rate_vs_sep(fill,rate_name1,"Raw","Nom",rate_name2,"Raw","Nom",2,1,iBC);    
  }
}

  /*
.L QA_comp_rate_vs_sep.C+
doFill("v0","t0",8379,20)
doFill("aodfddvtx","aodft0vtx",8379,20)
doFill("v0","aodfddvtx",8379,20)
doFill("t0","aodft0vtx",8379,20)

doFill("v0","t0",9128,20)
doFill("aodfddvtx","aodft0vtx",9128,20)
doFill("v0","aodfddvtx",9128,20)
doFill("t0","aodft0vtx",9128,20)

doFill("v0","t0",9644,20)
doFill("aodfddvtx","aodft0vtx",9644,20)
doFill("v0","aodfddvtx",9644,20)
doFill("t0","aodft0vtx",9644,20)

doFill("v0","t0",10298,20)
doFill("aodfddvtx","aodft0vtx",10298,20)
doFill("v0","aodfddvtx",10298,20)
doFill("t0","aodft0vtx",10298,20)

doFill("v0","t0",10824,20)
doFill("aodfddvtx","aodft0vtx",10824,20)
doFill("v0","aodfddvtx",10824,20)
doFill("t0","aodft0vtx",10824,20)

doFill("v0","t0",10782,32)
doFill("aodzn","aodft0vtx",10782,32)
doFill("v0","aodzn",10782,32)
doFill("t0","aodft0vtx",10782,32)
doFill("aodfddvtx","aodft0vtx",10782,32)

doFill("v0","t0",10802,40)
doFill("aodzn","aodft0vtx",10802,40)
doFill("aodzn","aodft0ceft0vtx",10802,40)
doFill("v0","aodzn",10802,40)
doFill("t0","aodft0vtx",10802,40)
doFill("t0","aodft0ceft0vtx",10802,40)
doFill("aodfddvtx","aodft0vtx",10802,40)

to merge use
cd aodPlots
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/v0_t0_8379_merged.pdf rate1_v0_rate2_t0_fill_8379*pdf
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/aodfddvtx_aodft0vtx_8379_merged.pdf rate1_aodfddvtx_rate2_aodft0vtx_fill_8379*pdf
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/v0_aodfddvtx_8379_merged.pdf rate1_v0_rate2_aodfddvtx_fill_8379*pdf
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/t0_aodft0vtx_8379_merged.pdf rate1_t0_rate2_aodft0vtx_fill_8379*pdf

gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/v0_t0_9128_merged.pdf rate1_v0_rate2_t0_fill_9128*pdf
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/aodfddvtx_aodft0vtx_9128_merged.pdf rate1_aodfddvtx_rate2_aodft0vtx_fill_9128*pdf
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/v0_aodfddvtx_9128_merged.pdf rate1_v0_rate2_aodfddvtx_fill_9128*pdf
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/t0_aodft0vtx_9128_merged.pdf rate1_t0_rate2_aodft0vtx_fill_9128*pdf

gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/v0_t0_9644_merged.pdf rate1_v0_rate2_t0_fill_9644*pdf
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/aodfddvtx_aodft0vtx_9644_merged.pdf rate1_aodfddvtx_rate2_aodft0vtx_fill_9644*pdf
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/v0_aodfddvtx_9644_merged.pdf rate1_v0_rate2_aodfddvtx_fill_9644*pdf
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/t0_aodft0vtx_9644_merged.pdf rate1_t0_rate2_aodft0vtx_fill_9644*pdf

gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/v0_t0_10298_merged.pdf rate1_v0_rate2_t0_fill_10298*pdf
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/aodfddvtx_aodft0vtx_10298_merged.pdf rate1_aodfddvtx_rate2_aodft0vtx_fill_10298*pdf
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/v0_aodfddvtx_10298_merged.pdf rate1_v0_rate2_aodfddvtx_fill_10298*pdf
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/t0_aodft0vtx_10298_merged.pdf rate1_t0_rate2_aodft0vtx_fill_10298*pdf

gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/v0_t0_10824_merged.pdf rate1_v0_rate2_t0_fill_10824*pdf
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/aodfddvtx_aodft0vtx_10824_merged.pdf rate1_aodfddvtx_rate2_aodft0vtx_fill_10824*pdf
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/v0_aodfddvtx_10824_merged.pdf rate1_v0_rate2_aodfddvtx_fill_10824*pdf
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/t0_aodft0vtx_10824_merged.pdf rate1_t0_rate2_aodft0vtx_fill_10824*pdf

gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/v0_t0_10782_merged.pdf rate1_v0_rate2_t0_fill_10782*pdf
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/aodzn_aodft0vtx_10782_merged.pdf rate1_aodzn_rate2_aodft0vtx_fill_10782*pdf
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/v0_aodzn_10782_merged.pdf rate1_v0_rate2_aodzn_fill_10782*pdf
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/t0_aodft0vtx_10782_merged.pdf rate1_t0_rate2_aodft0vtx_fill_10782*pdf
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/aodfddvtx_aodft0vtx_10782_merged.pdf rate1_aodfddvtx_rate2_aodft0vtx_fill_10782*pdf

gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/v0_t0_10802_merged.pdf rate1_v0_rate2_t0_fill_10802*pdf
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/aodzn_aodft0vtx_10802_merged.pdf rate1_aodzn_rate2_aodft0vtx_fill_10802*pdf
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/aodzna_odft0ceft0vtx_10802_merged.pdf rate1_aodzn_rate2_aodft0ceft0vtx_fill_10802*pdf
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/v0_aodzn_10802_merged.pdf rate1_v0_rate2_aodzn_fill_10802*pdf
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/t0_aodft0vtx_10802_merged.pdf rate1_t0_rate2_aodft0vtx_fill_10802*pdf
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/t0_aodft0ceft0vtx_10802_merged.pdf rate1_t0_rate2_aodft0ceft0vtx_fill_10802*pdf
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=../compSummary/aodfddvtx_aodft0vtx_10802_merged.pdf rate1_aodfddvtx_rate2_aodft0vtx_fill_10802*pdf
  */


