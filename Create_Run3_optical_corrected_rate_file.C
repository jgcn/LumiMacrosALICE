#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// Create a corrected rate file for a given
// trigger, using Ivan's files with the optical correction
//  ** WARNING: Directory structure and naming convention may change from fill to fill ...
//-------------------------------------------------------

void Create_Run3_optical_corrected_rate_file(Int_t Fill, const char* sys_opt)
// sys_opt = "" default, other options are for systematic studies:
// "+Q", "-Q", "+beta", "-beta", "dQip"
{

  // get name of files and set pointers to trees
  Set_input_file_names(Fill);
  Set_pointers_to_input_files_and_trees();

  // -- get number of separations for first scan
  //    (it is the same for both scans and for x and y)
  Int_t n_sep = FindNumberSeparations(1, 0);

  // get number of BCs 
  Int_t nIBC = GetNumberInteractingBunchCrossings();
  
  // to store names for rate files
  char *file_name_rate_x = new char[kg_string_size];  // input rates
  char *file_name_new_rate_x = new char[kg_string_size]; // new rates
  char *file_name_rate_y = new char[kg_string_size];  // input rates
  char *file_name_new_rate_y = new char[kg_string_size]; // new rates

  // create  file name with optical correction
  char *optical_file_name = new char[kg_string_size];
  sprintf(optical_file_name,"../Fill-%d/Corr-%d-sys/ODC-sys/ROOT/OPT-bbroot_%d_T0%s.root",
	  g_vdm_Fill,g_vdm_Fill,g_vdm_Fill,sys_opt);

  // open  file and get the trees
  TFile *optical_file = new TFile(optical_file_name);
  TTree *optical_tree = (TTree *) optical_file->Get("optical");
  optical_tree->ResetBranchAddresses();
  Double_t *optLL0 = new Double_t [n_sep*2];
  optical_tree->SetBranchAddress("optLL0",optLL0);
  Double_t LL000x;
  optical_tree->SetBranchAddress("LL000x",&LL000x);
  Double_t LL000y;
  optical_tree->SetBranchAddress("LL000y",&LL000y);

  // info in Ivan's trees seems to be organised
  // in an outer loop over scans, then over bc
  // and for each bc, half the info is x the other y sep.
  for(Int_t scan=0; scan<2;scan++) {
    cout << " doing  scan " << scan << endl;
    // prepare file names for the output and input

    // -- names of rate files
    sprintf(file_name_rate_x,"../Fill-%d/IntensityCorrFBCTRate_aodT0_x_Scan_%d.root",
	     g_vdm_Fill,scan);
    sprintf(file_name_rate_y,"../Fill-%d/IntensityCorrFBCTRate_aodT0_y_Scan_%d.root",
	     g_vdm_Fill,scan);
    // -- names of new rate files
    sprintf(file_name_new_rate_x,"../Fill-%d/OpticalIntensityCorrFBCTRate_aodT0_x_Scan_%d.root",
	     g_vdm_Fill,scan);
    sprintf(file_name_new_rate_y,"../Fill-%d/OpticalIntensityCorrFBCTRate_aodT0_y_Scan_%d.root",
	     g_vdm_Fill,scan);


    // open file with rates in x
    TFile *RateFile_x = new TFile(file_name_rate_x);
    TTree *rate_tree_x = (TTree *) RateFile_x->Get("Rate");
    Double_t *rate_x = new Double_t[n_sep];
    Double_t *rate_error_x = new Double_t[n_sep];  
    rate_tree_x->SetBranchAddress("rate",rate_x);
    rate_tree_x->SetBranchAddress("rate_error",rate_error_x);   
    // open file with rates in y
    TFile *RateFile_y = new TFile(file_name_rate_y);
    TTree *rate_tree_y = (TTree *) RateFile_y->Get("Rate");
    Double_t *rate_y = new Double_t[n_sep];
    Double_t *rate_error_y = new Double_t[n_sep];  
    rate_tree_y->SetBranchAddress("rate",rate_y);
    rate_tree_y->SetBranchAddress("rate_error",rate_error_y);
    
    // create trees with new rates
    TFile *RateFile_new_x = new TFile(file_name_new_rate_x,"recreate");
    TTree *new_rate_x_tree = new TTree("Rate","Rate");
    Double_t *new_rate_x = new Double_t[n_sep];
    Double_t *new_rate_error_x = new Double_t[n_sep];
    char *txt_tmp = new char[kg_string_size];
    sprintf(txt_tmp,"rate[%d]/D",n_sep);
    new_rate_x_tree->Branch("rate",new_rate_x,txt_tmp);
    sprintf(txt_tmp,"rate_error[%d]/D",n_sep);
    new_rate_x_tree->Branch("rate_error",new_rate_error_x,txt_tmp);
    TFile *RateFile_new_y = new TFile(file_name_new_rate_y,"recreate");
    TTree *new_rate_y_tree = new TTree("Rate","Rate");
    Double_t *new_rate_y = new Double_t[n_sep];
    Double_t *new_rate_error_y = new Double_t[n_sep];
    sprintf(txt_tmp,"rate[%d]/D",n_sep);
    new_rate_y_tree->Branch("rate",new_rate_y,txt_tmp);
    sprintf(txt_tmp,"rate_error[%d]/D",n_sep);
    new_rate_y_tree->Branch("rate_error",new_rate_error_y,txt_tmp);

    // loop over bunches
    for(Int_t k=0;k<nIBC;k++) {
      optical_tree->GetEntry(k+nIBC*scan);
      rate_tree_x->GetEntry(k);
      rate_tree_y->GetEntry(k);      
      for(Int_t isep = 0;isep<n_sep;isep++) {
	new_rate_x[isep] = LL000x*rate_x[isep]/optLL0[isep];
	new_rate_y[isep] = LL000y*rate_y[isep]/optLL0[isep+n_sep];	
	new_rate_error_x[isep] = LL000x*rate_error_x[isep]/optLL0[isep];
	new_rate_error_y[isep] = LL000y*rate_error_y[isep]/optLL0[isep+n_sep];
      }
      RateFile_new_x->cd();
      new_rate_x_tree->Fill();
      RateFile_new_y->cd();
      new_rate_y_tree->Fill();
    } // loop over BCid
    // save
    RateFile_new_x->cd();
    new_rate_x_tree->Write();
    RateFile_new_x->Close();
    RateFile_new_y->cd();
    new_rate_y_tree->Write();
    RateFile_new_y->Close();

    // clean
    delete [] rate_x;
    delete [] rate_error_x;    
    delete [] rate_y;
    delete [] rate_error_y;    
    delete [] new_rate_x;
    delete [] new_rate_error_x;    
    delete [] new_rate_y;
    delete [] new_rate_error_y;    
  } // end loop over scans
  
  // clean up
  delete [] file_name_rate_x;
  delete [] file_name_new_rate_x;
  delete [] file_name_rate_y;
  delete [] file_name_new_rate_y;
  delete [] optical_file_name;
  delete [] optLL0;

}
