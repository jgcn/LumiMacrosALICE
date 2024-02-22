#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// Macro to print the start, end and duration of the info in Ivan's trees
//-------------------------------------------------------

void QA_print_time_start_end(Int_t Fill)
{
	// initialize
	Set_input_file_names(Fill);
	Set_pointers_to_input_files_and_trees();

	// get the histo
	TH1D* h = (TH1D*) g_vdm_File->Get("tstat");

	// get the times
	int64_t startTime = h->GetBinContent(1);
	int64_t endTime = h->GetBinContent(2);

	// print out
	cout << " Start = " << startTime << ", end = " << endTime
	     <<", duration = " << (endTime-startTime) << endl;
}
