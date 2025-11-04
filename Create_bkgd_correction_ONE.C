#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// create a file containing background-correction factors equal to one
// that is, no background correction
//-------------------------------------------------------


//-------------------------------------------------------
void createOneFile(int scan, int direction)
// direction 0 = x, 1 = y
{
  // create output file	
  const char *directionName[2] = {"x","y"};		
  const char* file_name = Form("../Fill-%i/BkgdCorr_ONE_%s_Scan_%i.root",
			       g_vdm_Fill,directionName[direction], scan);
  TFile *CorrFile = new TFile(file_name,"recreate");

  // create output tree
  Int_t n_separations = FindNumberSeparations(direction+1, scan);
  Double_t *correction = new Double_t[n_separations];
  Double_t *correction_error = new Double_t[n_separations];  
  TTree *correction_tree = new TTree("BkgdCorr","BkgdCorr");
  correction_tree->Branch("correction",correction,Form("correction[%d]/D",n_separations));
  correction_tree->Branch("correction_error",correction_error,Form("correction_error[%d]/D",n_separations));
  
  // fill the file
  Int_t nBC = GetNumberInteractingBunchCrossings();
  std::cout << " bcs " << nBC << " separations " << n_separations << std::endl;
  for(int iBC=0;iBC<nBC;iBC++) {
    for(int i=0;i<n_separations;i++) {
      // initialise
      correction_error[i]=0;
      correction[i]=1;
    } // loop over separations
    correction_tree->Fill();
  } // loop over BCs

  // save tree
  CorrFile->cd();
  correction_tree->SetDirectory(CorrFile);
  correction_tree->Write();
  CorrFile->Close();

}

//-------------------------------------------------------
void Create_bkgd_correction_ONE(Int_t Fill) 
{
  // get name of files and set pointers to trees
  Set_input_file_names(Fill);
  Set_pointers_to_input_files_and_trees();
  
  // create files for all scans
  for (Int_t i=0;i<g_n_Scans_in_Fill;i++)
    {
      createOneFile(i,0);
      createOneFile(i,1);      
    }
  
  return;
}
