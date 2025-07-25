#ifndef __INPUT_FROM_USER_HH__
#define __INPUT_FROM_USER_HH__

//-------------------------------------------------------
// here the user defines the variables
// to be used as input for the analysis
//-------------------------------------------------------

void Set_input_file_names(Int_t Fill)
{
  //  	if (Fill == 8379) // pp @ 13.6 TeV Nov 10, 2022
  //8379-> pp @ 13.6 TeV Nov 10 2022, 9128-> pp @ 13.6 TeV Sep 7-8 2023, 9644-> pp @ 13.6 TeV May 18 2024
  // 10298-> pp @ 5.36 TeV, Oct-28-24
  // 10782 ->pO , 10802 -> OO, 10813 -> NeNe
  // 10813 -> pp 2025
       	if (Fill == 8379 || Fill ==9128 || Fill ==9644 || Fill ==10298
	    || Fill ==10782 || Fill ==10802 || Fill ==10813 || Fill ==10824)  
	{
		// set fill and number of scans in the fill
		g_vdm_Fill = Fill;
		g_n_Scans_in_Fill = 2;

		const char* CUT = "DUMMY";
		//cout <<" Following V0/T0 timing cut is being used: " <<CUT <<endl;
		// set name of input files
		sprintf(g_Input_vdm_File,      "../Fill-%d/vdm_time_%d_%s_1_v3.root",      g_vdm_Fill, g_vdm_Fill, CUT);
		sprintf(g_Input_vdm_BPTX_File, "../Fill-%d/vdm_time_%d_%s_1_v3-BPTX.root", g_vdm_Fill, g_vdm_Fill, CUT);

		// charge of beams
		gBeamA = 1; // proton
		gBeamB = 1; // proton
		if (Fill == 10782) {gBeamB = 8;}
		if (Fill == 10802) {gBeamA = gBeamB = 8;}
		if (Fill == 10813) {gBeamA = gBeamB = 10;g_n_Scans_in_Fill = 1;}
				
	}
	else
	{
		cout << " Fill " << Fill << " not know. Bye " << endl;
		exit(-100);
	}

	return;
}
#endif
