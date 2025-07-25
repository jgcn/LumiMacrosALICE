#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"
#include "FitUtils.h"
#include "Plotting.h"

//-------------------------------------------------------
// Macro to produce quality-assurance (QA) plots related to
// the cross section
//-------------------------------------------------------

TGraphErrors* GetGraphHxHy(const char *rate_name, const char *rate_type,
                  const char *sep_type, const char *intensity_type, Int_t fit_type, Int_t scan)
{
    
	// get number of bunch crossings
	Int_t nIBC = GetNumberInteractingBunchCrossings();
	Int_t Bunches[nIBC]; //kimc
	GetBunchIndices(Bunches);
	
	//Get bad bunches list, kimc
	SetBCBlacklists();

	// --------------------
	// get intensity file and tree
	//---------------------
	// open the correct file
	char filename[120];
	sprintf(filename,"../Fill-%d/%s_Scan_%d.root",g_vdm_Fill,intensity_type,scan);
	TFile *IntensityFile = new TFile(filename);

	// get the info from the tree
	Double_t *bunch_intensity_1 = new Double_t[nIBC];
	Double_t *bunch_intensity_2 = new Double_t[nIBC];
	Double_t cf_dcct_1;
	Double_t cf_dcct_2;

	TTree *intensity_tree = (TTree *) IntensityFile->Get("Bunch_Intensity");
	intensity_tree->ResetBranchAddresses();
	intensity_tree->SetBranchAddress("cf_dcct_1",&cf_dcct_1);
	intensity_tree->SetBranchAddress("cf_dcct_2",&cf_dcct_2);    
	intensity_tree->SetBranchAddress("bunch_intensity_1",bunch_intensity_1);
	intensity_tree->SetBranchAddress("bunch_intensity_2",bunch_intensity_2);
	intensity_tree->GetEntry(0);

        // reserve space
	Double_t *N1N2_all = new Double_t [nIBC];
     
    // first get the files and trees
    // --> create hx/hy file names
    char *hx_file_name = new char[kg_string_size];
    sprintf(hx_file_name,"../Fill-%d/hxy_%sRate_%s_%sSep_x_Scan_%d_Fit_%s.root",
            g_vdm_Fill,rate_type,rate_name,sep_type,scan,g_fit_model_name[fit_type]);
    char *hy_file_name = new char[kg_string_size];
    sprintf(hy_file_name,"../Fill-%d/hxy_%sRate_%s_%sSep_y_Scan_%d_Fit_%s.root",
            g_vdm_Fill,rate_type,rate_name,sep_type,scan,g_fit_model_name[fit_type]);

    // --> open files
    TFile *hx_file = new TFile(hx_file_name);
    TFile *hy_file = new TFile(hy_file_name);

    // --> get the trees
    TTree *hx_tree = (TTree *) hx_file->Get("AreaRate");
    TTree *hy_tree = (TTree *) hy_file->Get("AreaRate");
   
    // Next step: prepare variables to store the info for each bunch crossing
    // --> variables
    Double_t chi2_dof_x;
    Double_t chi2_dof_y;
    Double_t *area_x = new Double_t[2]; // area and its error
    Double_t *rate_zero_x = new Double_t[2]; // rate at zero and its error
    Double_t *area_y = new Double_t[2]; // area and its error
    Double_t *rate_zero_y = new Double_t[2]; // rate at zero and its error
    // --> set branches for hx, hy
    hx_tree->ResetBranchAddresses();
    hy_tree->ResetBranchAddresses();
    hx_tree->SetBranchAddress("chi2_dof",&chi2_dof_x);
    hy_tree->SetBranchAddress("chi2_dof",&chi2_dof_y);
    hx_tree->SetBranchAddress("area",area_x);
    hy_tree->SetBranchAddress("area",area_y);
    hx_tree->SetBranchAddress("rate_zero",rate_zero_x);
    hy_tree->SetBranchAddress("rate_zero",rate_zero_y);
      
    // define graph    
	TGraphErrors* gr = new TGraphErrors();
	//gr->SetMarkerColor(0);
	
	// fill info
	Int_t n = 0;
    
    // fill info
    for(Int_t i=0;i<nIBC;i++) {
        //Get info
        hx_tree->GetEntry(i);
        hy_tree->GetEntry(i);
        
        if (OnBCBlacklist(g_vdm_Fill, Bunches[i])) //Rule out bad bunches, kimc
		{
			cout <<Form("Bad bunch detected, rule it out: %i (index %i)\n", Bunches[i], i);
			continue;
		}
	
         Double_t N1N2_all = cf_dcct_1 * bunch_intensity_1[i] * cf_dcct_2 * bunch_intensity_2[i];
               N1N2_all = N1N2_all / (1.E18);   //July 1, scale down x axis by 10^18 for PN not plot - kimc
 
   //   if (!(UseBunchCrossing( i))) continue;
     // if (chi2_dof_x <= 0 || chi2_dof_y <= 0) continue;
     
        if (!(chi2_dof_x > 0 && chi2_dof_y > 0)) continue;
            //fill histograms
            Double_t hxhy = GetHxHy(area_x[0],area_y[0],rate_zero_x[0],rate_zero_y[0]);
            Double_t hxhy_err = GetHxHyerr(area_x[0],area_x[1],area_y[0],area_y[1],rate_zero_x[0],rate_zero_x[1],rate_zero_y[0],rate_zero_y[1]);
           
     //       cout <<"bin: "<< i<< " " << hxhy<< " +/-" <<  hxhy_err << endl;
		
	      gr->SetPoint(n, N1N2_all, hxhy);
              gr->SetPointError(n, 0, hxhy_err); // sin error en X, solo en Y
             cout <<"Point: "<< n << "  N1N2 = " << N1N2_all << "  hxhy = " << hxhy << endl;
                n++;               
    }

    return gr;
    
}

void QA_hxhy_V0T0_N1N2(Int_t Fill, const char *rate_name1, const char *rate_name2, const char *rate_type, const char *sep_type,
              const char *intensity_type, Int_t fit_type,
              Bool_t save = kTRUE)
{
    // initialize
    Set_input_file_names(Fill);
    Set_pointers_to_input_files_and_trees();
    TString name = Form("Ratio of h_{x}h_{y} #frac{%s}{%s}",rate_name1,rate_name2);
    
	//Plot for public note, June 25
    TGraphErrors* gr_T0_0 = GetGraphHxHy(rate_name1, rate_type, sep_type, intensity_type, fit_type, 0);
    TGraphErrors* gr_T0_1 = GetGraphHxHy(rate_name1, rate_type, sep_type, intensity_type, fit_type, 1);
    TGraphErrors* gr_V0_0 = GetGraphHxHy(rate_name2, rate_type, sep_type, intensity_type, fit_type, 0);
    TGraphErrors* gr_V0_1 = GetGraphHxHy(rate_name2, rate_type, sep_type, intensity_type, fit_type, 1);

    // Crear ratio = T0/V0 punto a punto
auto DivideGraphsWithErrors = [](TGraphErrors* num, TGraphErrors* den) -> TGraphErrors* {
    TGraphErrors* gr_ratio = new TGraphErrors();
    for (int i = 0; i < num->GetN(); i++) {
        double x1, y1, ey1;
        num->GetPoint(i, x1, y1);
        ey1 = num->GetErrorY(i);
        bool found = false;
        for (int j = 0; j < den->GetN(); j++) {
            double x2, y2, ey2;
            den->GetPoint(j, x2, y2);
            ey2 = den->GetErrorY(j);
            if (fabs(x1 - x2) < 1e-6 && y2 != 0) {
                double ratio = y1 / y2;
                double eratio = ratio * sqrt(pow(ey1 / y1, 2) + pow(ey2 / y2, 2));
                gr_ratio->SetPoint(gr_ratio->GetN(), x1, ratio);
                gr_ratio->SetPointError(gr_ratio->GetN()-1, 0, eratio);
                found = true;
                break;
            }
        }
        if (!found)
            std::cerr << "Warning: No match in denominator for x = " << x1 << std::endl;
    }
    return gr_ratio;
};
    TGraphErrors* ratio[2];
    ratio[0] = DivideGraphsWithErrors(gr_T0_0, gr_V0_0);
    ratio[1] = DivideGraphsWithErrors(gr_T0_1, gr_V0_1);
    
   // Plot
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    TCanvas* c1 = new TCanvas("c1_RatioT0V0", "", 1600, 1000);
    c1->Divide(1, 2);

    for (int a = 0; a < 2; a++) {
        c1->cd(a + 1);
        ratio[a]->SetTitle(Form("Scan %d;N1N2 [10^{18}];%s", a, name.Data()));
        ratio[a]->SetMarkerStyle(20);
        ratio[a]->SetMarkerSize(1.4);
        ratio[a]->GetYaxis()->SetRangeUser(0.9, 1.04);
        ratio[a]->Draw("AP");

        TF1* f1 = new TF1(Form("F1_ratio%d", a), "pol0", 0, 9999);
        f1->SetLineStyle(2);
        ratio[a]->Fit(f1, "Q");

     //   TLegend *L = new TLegend(0.15, 0.7, 0.45, 0.85);
     //   L->AddEntry((TObject*)0, Form("Mean = %.3f #pm %.3f", f1->GetParameter(0), f1->GetParError(0)), "");
     //   L->Draw();
        
        	TLegend *L1 = new TLegend(0.15, 0.7, 0.35, 0.825);
		L1->SetMargin(0);
		L1->SetBorderSize(0);
		L1->SetFillStyle(3001);
		L1->AddEntry((TObject*)0, "ALICE", "");
		L1->AddEntry((TObject*)0, "pp #sqrt{s} = 13.6 TeV", "");
		L1->Draw();

/*
		TLegend *L2 = new TLegend(0.15, 0.6, 0.35, 0.7);
		L2->SetMargin(0);
		L2->SetBorderSize(0);
		L2->SetFillStyle(3001);
		L2->AddEntry((TObject*)0, Form("Scan %i", a), ""); 
		L2->Draw();
*/
		TLegend *L3 = new TLegend(0.675, 0.72, 0.875, 0.85);
		L3->SetMargin(0.1);
		L3->SetNColumns(2);
		L3->AddEntry((TObject*)0, "#chi^{2}/NDF", "");
		L3->AddEntry((TObject*)0, Form("%5.2f/%i", f1->GetChisquare(), f1->GetNDF()), "");
		L3->AddEntry((TObject*)0, "Mean", "");
		L3->AddEntry((TObject*)0, Form("%5.3f #pm %4.3f", f1->GetParameter(0), f1->GetParError(0)), "");
		L3->Draw();
        
        
    }
    if (save)
        c1->Print(Form("c5a_Fill%i_Ratiohxhy_V0T0_N1N2.eps", Fill));
}





