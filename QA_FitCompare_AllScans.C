#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"
#include "TPaveText.h"
#include "TMultiGraph.h"
#include "FitUtils.h"

void QA_FitCompare_AllScans(Int_t Fill, const char *rate_name, const char *rate_type,
                            const char *sep_type) {
    gROOT->SetBatch(kTRUE);  // no mostrar gráficos en pantalla
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    const char* fit_names[] = {"GP2", "GP6", "G", "DG"};
    int fit_types[] = {0, 1, 2, 4};
    Color_t colors[] = {kRed, kBlue, kGreen+2, kMagenta};

    Set_input_file_names(Fill);
    Set_pointers_to_input_files_and_trees();
    

    // PDF multipágina único
    TString pdfFile = Form("FitCompare_AllBCs_Fill%d_%s.pdf", Fill, rate_name);
    TString pdfOpen = pdfFile + "[";
    TString pdfClose = pdfFile + "]";
    TCanvas *c = new TCanvas("c", "Rate vs Separation", 800, 600);
    c->Print(pdfOpen);

    // Recorremos todos los casos: X0, X1, Y0, Y1
    for (int scan_type = 1; scan_type <= 2; ++scan_type) {
        for (int scan = 0; scan <= 1; ++scan) {

            // Archivos de entrada
            char rate_file_name[300], sep_file_name[300];
            if (scan_type == 1) {
                sprintf(rate_file_name, "../Fill-%d/%sRate_%s_x_Scan_%d.root", Fill, rate_type, rate_name, scan);
                sprintf(sep_file_name, "../Fill-%d/%sSep_x_Scan_%d.root", Fill, sep_type, scan);
            } else {
                sprintf(rate_file_name, "../Fill-%d/%sRate_%s_y_Scan_%d.root", Fill, rate_type, rate_name, scan);
                sprintf(sep_file_name, "../Fill-%d/%sSep_y_Scan_%d.root", Fill, sep_type, scan);
            }

            TFile *rate_file = TFile::Open(rate_file_name);
            TFile *sep_file = TFile::Open(sep_file_name);
            if (!rate_file || !sep_file || rate_file->IsZombie() || sep_file->IsZombie()) {
                std::cerr << "Error: no se pudieron abrir " << rate_file_name 
                          << " o " << sep_file_name << std::endl;
                continue;
            }

            TTree *rate_tree = (TTree*) rate_file->Get("Rate");
            TTree *sep_tree = (TTree*) sep_file->Get("Separations");
            if (!rate_tree || !sep_tree) {
                std::cerr << "Error: no se encontraron los árboles." << std::endl;
                continue;
            }

            Int_t nBCs = rate_tree->GetEntries();
            Int_t n_sep = FindNumberSeparations(scan_type, scan);

            Double_t *rate = new Double_t[n_sep];
            Double_t *rate_error = new Double_t[n_sep];
            Double_t *sep = new Double_t[n_sep];

            rate_tree->SetBranchAddress("rate", rate);
            rate_tree->SetBranchAddress("rate_error", rate_error);
            sep_tree->SetBranchAddress("separation", sep);

            // Loop sobre todos los BCs
            for (int bc = 0; bc < nBCs; ++bc) {
                rate_tree->GetEntry(bc);
                sep_tree->GetEntry(bc);

                c->Clear();
                TMultiGraph *mg = new TMultiGraph();
                TLegend *leg = new TLegend(0.65, 0.65, 0.88, 0.88);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextSize(0.03);

                // Datos
                TGraphErrors *grData = new TGraphErrors(n_sep, sep, rate, nullptr, rate_error);
                grData->SetMarkerStyle(20);
                grData->SetMarkerColor(kBlack);
                grData->SetLineColor(kBlack);
                mg->Add(grData, "P");
                leg->AddEntry(grData, "Fit method", "p");

                // Ajustes
                for (int i = 0; i < 4; ++i) {
                    Int_t fit_type = fit_types[i];
                    const char* fit_label = fit_names[i];

                    Double_t area[2], rate_zero[2], chi2_dof;
                    Int_t npar = Get_number_par(fit_type);
                    Double_t *par = new Double_t[npar];
                    Double_t *par_err = new Double_t[npar];

                    TString cname = Form("Fit_%s_Fill%d_Scan%d_BC%d_Rate%s", fit_label, Fill, scan, bc, rate_name);
		            const char* corrName = Form("CorrelationMatrix_i%i", bc);
                    TH2D *hCorr = new TH2D(corrName, corrName, npar, 0, npar, npar, 0, npar); //can be written to file, if needed
                    chi2_dof = Fit_rate_separation(n_sep, sep, rate, rate_error, fit_type,
                                                   area, rate_zero, par, par_err, hCorr, cname.Data());

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
                    mg->Add(gr, "L");

                    if (chi2_dof >= 0)
                        leg->AddEntry(gr, Form("%s (#chi^{2}/ndf=%.2f)", fit_label, chi2_dof), "l");
                    else
                        leg->AddEntry(gr, Form("%s (fit failed)", fit_label), "l");

                    delete[] par;
                    delete[] par_err;
                    delete[] fit_vals;
                }



                // Título con etiquetas
                mg->SetTitle(";Separation (mm);Rate [Hz]");  // solo ejes
                c->cd()->SetLogy();
                mg->Draw("ALP");
                leg->Draw();

                // Titulo visible
                TPaveText *pt = new TPaveText(0.15, 0.92, 0.85, 0.99, "NDC");
                pt->SetFillColor(0);
                pt->SetBorderSize(0);
                pt->SetTextAlign(22);
                pt->SetTextFont(42);
                pt->SetTextSize(0.04);
                pt->AddText(Form("Fill %d |  Rate %s  |  Scan %s%d  |  BC %d ", 
                 Fill, rate_name, (scan_type==1?"X":"Y"), scan, bc));
                pt->Draw();

                c->Print(pdfFile); // agregar hoja al pdf
                
            }

            delete[] rate;
            delete[] rate_error;
            delete[] sep;
        }
    }

    // cerrar pdf
    c->Print(pdfClose);
    std::cout << "PDF generado: " << pdfFile << std::endl;
}

