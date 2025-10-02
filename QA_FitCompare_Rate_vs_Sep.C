
#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"
#include "TPaveText.h"
#include "TMultiGraph.h"
#include "FitUtils.h"

void QA_FitCompare_Rate_vs_Sep(Int_t Fill, const char *rate_name, const char *rate_type,
                            const char *sep_type, Int_t scan_type, Int_t scan, Int_t bc) {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

   // const char* fit_names[] = {"GP6", "GP2", "DG", "Num"};
   //  const char* fit_names[] = {"GP2", "GP6","G","NUM", "DG"};
          const char* fit_names[] = {"GP2", "GP6","G", "DG"};
    int fit_types[] = {0, 1, 2, 4};
    Color_t colors[] = {kRed, kBlue, kGreen+2, kMagenta};

    Set_input_file_names(Fill);
    Set_pointers_to_input_files_and_trees();

    char rate_file_name[300], sep_file_name[300];
    if (scan_type == 1) {
        sprintf(rate_file_name, "../Fill-%d/%sRate_%s_x_Scan_%d.root", Fill, rate_type, rate_name, scan);
        sprintf(sep_file_name, "../Fill-%d/%sSep_x_Scan_%d.root", Fill, sep_type, scan);
    } else {
        sprintf(rate_file_name, "../Fill-%d/%sRate_%s_y_Scan_%d.root", Fill, rate_type, rate_name, scan);
        sprintf(sep_file_name, "../Fill-%d/%sSep_y_Scan_%d.root", Fill, sep_type, scan);
    }

    TFile *rate_file = new TFile(rate_file_name);
    TFile *sep_file = new TFile(sep_file_name);
    if (!rate_file || !sep_file) {
        std::cerr << "Error open files.\n";
        return;
    }

    TTree *rate_tree = (TTree*) rate_file->Get("Rate");
    TTree *sep_tree = (TTree*) sep_file->Get("Separations");

    Int_t n_sep = FindNumberSeparations(scan_type, scan);
    Double_t *rate = new Double_t[n_sep];
    Double_t *rate_error = new Double_t[n_sep];
    Double_t *sep = new Double_t[n_sep];

    rate_tree->SetBranchAddress("rate", rate);
    rate_tree->SetBranchAddress("rate_error", rate_error);
    sep_tree->SetBranchAddress("separation", sep);

    rate_tree->GetEntry(bc);
    sep_tree->GetEntry(bc);

    TCanvas *c = new TCanvas("c", "Rate vs Separation - Fitted", 800, 600);
    TMultiGraph *mg = new TMultiGraph();
    TLegend *leg = new TLegend(0.65, 0.65, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.03);

    TGraphErrors *grData = new TGraphErrors(n_sep, sep, rate, nullptr, rate_error);
    grData->SetMarkerStyle(20);
    grData->SetMarkerColor(kBlack);
    grData->SetLineColor(kBlack);
    grData->SetTitle("Datos");
    mg->Add(grData, "P");
    leg->AddEntry(grData, "Datos", "p");

    for (int i = 0; i < 4; ++i) {
        Int_t fit_type = fit_types[i];
        const char* fit_label = fit_names[i];

        Double_t area[2], rate_zero[2], chi2_dof;
        Int_t npar = Get_number_par(fit_type);
        Double_t *par = new Double_t[npar];
        Double_t *par_err = new Double_t[npar];

        TString cname = Form("Fit_%s_Fill%d_Scan%d_BC%d", fit_label, Fill, scan, bc);
        chi2_dof = Fit_rate_separation(n_sep, sep, rate, rate_error, fit_type,
                                       area, rate_zero, par, par_err, cname.Data());

        Double_t *fit_vals = new Double_t[n_sep];
        for (int j = 0; j < n_sep; ++j) {
            Double_t x[1] = {sep[j]};
            if      (fit_type == 0) fit_vals[j] = fit_GP2(x, par);
            else if (fit_type == 1) fit_vals[j] = fit_GP6(x, par);
            else if (fit_type == 2) fit_vals[j] = fit_G(x, par);
            else if (fit_type == 4) fit_vals[j] = fit_DG(x, par);
            else                    fit_vals[j] = NAN;
        }

        TGraph *gr = new TGraph(n_sep, sep, fit_vals);
        gr->SetLineColor(colors[i]);
        gr->SetLineWidth(2);
        gr->SetTitle(fit_label);
        mg->Add(gr, "L");

        if (chi2_dof >= 0)
            leg->AddEntry(gr, Form("%s (#chi^{2}/ndf=%.2f)", fit_label, chi2_dof), "l");
        else
            leg->AddEntry(gr, Form("%s (fit failed)", fit_label), "l");

        delete[] par;
        delete[] par_err;
        delete[] fit_vals;
    }
     
     c->cd(); // important!
     mg->SetTitle(Form("Rate vs Separation (%s - BC %d);Separation (mm);Rate [Hz]", rate_name, bc));
     c->cd()->SetLogy();
    mg->Draw("ALP"); // A=axes, L=lines, P=points
    leg->Draw();
    c->Update();

    TString fname = Form("FitCompare_RateVsSep_Fill%d_ScanT%d_Scan%d_BC%d.pdf", Fill, scan_type, scan, bc);
    c->SaveAs(fname);

    delete[] rate;
    delete[] rate_error;
    delete[] sep;
}
