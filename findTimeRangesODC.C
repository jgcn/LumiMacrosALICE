#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

//-------------------------------------------------------
// Find the time ranges that can be used to fit orbit-drift data
// -- time variables are obtained by running QA_print_timestamps.C
//-------------------------------------------------------

// varibles to store the time ranges
// before/middle of/after the sacan
int64_t startTime[3];
int64_t endTime[3];

void initRanges(int fill, int scan) {
  int offset = 100; // seconds
  if (fill == 8379 && scan == 0) {
    startTime[0] = 1668080744-offset;
    endTime[0] = 1668080744;
    startTime[1] = 1668081675;
    endTime[1] = 1668081877;
    startTime[2] = 1668082810;
    endTime[2] = 1668082810+offset;
  } else   if (fill == 8379 && scan == 1) {
    startTime[0] = 1668083731-offset;
    endTime[0] = 1668083731;
    startTime[1] = 1668084664;
    endTime[1] = 1668084862;
    startTime[2] = 1668085795;
    endTime[2] = 1668085795+offset;
  }
  if (fill == 9128 && scan == 0) {
    startTime[0] = 1694118912-offset;
    endTime[0] = 1694118912;
    startTime[1] = 1694119841;
    endTime[1] = 1694120095;
    startTime[2] = 1694121024;
    endTime[2] = 1694121024+offset;
  } else   if (fill == 9128 && scan == 1) {
    startTime[0] = 1694121414-offset;
    endTime[0] = 1694121414;
    startTime[1] = 1694122343;
    endTime[1] = 1694122619;
    startTime[2] = 1694123548;
    endTime[2] = 1694123548+offset;
  }
  if (fill == 9240 && scan == 0) {
    startTime[0] = 1696965052-offset;
    endTime[0] = 1696965052;
    startTime[1] = 1696965884;
    endTime[1] = 1696966032;
    startTime[2] = 1696966865;
    endTime[2] = 1696966865+offset;
  } else   if (fill == 9240 && scan == 1) {
    startTime[0] = 1696967005-offset;
    endTime[0] = 1696967005;
    startTime[1] = 1696967838;
    endTime[1] = 1696967994;
    startTime[2] = 1696968825;
    endTime[2] = 1696968825+offset;
  }
  if (fill == 9644 && scan == 0) {
    startTime[0] = 1716074080-offset;
    endTime[0] = 1716074080;
    startTime[1] = 1716074983;
    endTime[1] = 1716075207;
    startTime[2] = 1716076108;
    endTime[2] = 1716076108+offset;
  } else   if (fill == 9644 && scan == 1) {
    startTime[0] = 1716076416-offset;
    endTime[0] = 1716076416;
    startTime[1] = 1716077319;
    endTime[1] = 1716077559;
    startTime[2] = 1716078460;
    endTime[2] = 1716078460+offset;
  }
  if (fill == 10298 && scan == 0) {
    startTime[0] = 1730149754-offset;
    endTime[0] = 1730149754;
    startTime[1] = 1730150610;
    endTime[1] = 1730150717;
    startTime[2] = 1730151571;
    endTime[2] = 1730151571+offset;
  } else   if (fill == 10298 && scan == 1) {
    startTime[0] = 1730151710-offset;
    endTime[0] = 1730151710;
    startTime[1] = 1730152566;
    endTime[1] = 1730152673;
    startTime[2] = 1730153529;
    endTime[2] = 1730153529+offset;
  } 
  if (fill == 10782 && scan == 0) {
    startTime[0] = 1751423005-offset;
    endTime[0] = 1751423005;
    startTime[1] = 1751423700;
    endTime[1] = 1751423832;
    startTime[2] = 1751424529;
    endTime[2] = 1751424529+offset;
  } else   if (fill == 10782 && scan == 1) {
    startTime[0] = 1751426469-offset;
    endTime[0] = 1751426469;
    startTime[1] = 1751427166;
    endTime[1] = 1751427288;
    startTime[2] = 1751427985;
    endTime[2] = 1751427985+offset;
  } 
  if (fill == 10802 && scan == 0) {
    startTime[0] = 1751688982-offset;
    endTime[0] = 1751688982;
    startTime[1] = 1751689670;
    endTime[1] = 1751689778;
    startTime[2] = 1751690465;
    endTime[2] = 1751690465+offset;
  } else   if (fill == 10802 && scan == 1) {
    startTime[0] = 1751690567-offset;
    endTime[0] = 1751690567;
    startTime[1] = 1751691252;
    endTime[1] = 1751691344;
    startTime[2] = 1751692031;
    endTime[2] = 1751692031+offset;
  } 
  
}

void findTimeRangesODC(int fill)
{
  // get name of files and set pointers to trees
  Set_input_file_names(fill);
  Set_pointers_to_input_files_and_trees();
  
  // create files for all scans
  for(Int_t scan=0;scan<g_n_Scans_in_Fill;scan++)
    {
      initRanges(fill, scan);

      // set up input tree
      Int_t aqflag;
      Double_t time;
      g_vdm_Tree->ResetBranchAddresses();
      g_vdm_Tree->SetBranchAddress("aqflag",&aqflag);
      g_vdm_Tree->SetBranchAddress("time",&time);

      // store found times
      int64_t startTimeAQ[3];
      int64_t endTimeAQ[3];
      for(int i=0;i<3;i++) startTimeAQ[i]=endTimeAQ[i]=-1;
      // loop over the tree
      int nEntries = g_vdm_Tree->GetEntries();
      for (int i=0; i<nEntries;i++) {
	g_vdm_Tree->GetEntry(i);
	if (aqflag==0) continue;
	uint64_t currentTime = (uint64_t) time;
	if (startTimeAQ[0] == -1 && startTime[0]<currentTime) startTimeAQ[0] = currentTime;
	if (startTime[0]<currentTime && endTime[0]>currentTime) endTimeAQ[0] =  currentTime;
	if (startTimeAQ[1] == -1 && startTime[1]<currentTime) startTimeAQ[1] = currentTime;
	if (startTime[1]<currentTime && endTime[1]>currentTime) endTimeAQ[1] =  currentTime;
	if (startTimeAQ[2] == -1 && startTime[2]<currentTime) startTimeAQ[2] = currentTime;
	if (startTime[2]<currentTime && endTime[2]>currentTime) endTimeAQ[2] =  currentTime;
	if (endTime[2]<currentTime) break;
      }

      // print times
      std::cout << " Fill " << fill << ", scan " << scan << std::endl;
      std::cout << "   startTimeAQ[0] = " << startTimeAQ[0] << ";" << std::endl;
      std::cout << "   startTimeAQ[1] = " << startTimeAQ[1] << ";" << std::endl;
      std::cout << "   startTimeAQ[2] = " << startTimeAQ[2] << ";" << std::endl;
      std::cout << "   endTimeAQ[0] = " << endTimeAQ[0] << ";" << std::endl;
      std::cout << "   endTimeAQ[1] = " << endTimeAQ[1] << ";" << std::endl;
      std::cout << "   endTimeAQ[2] = " << endTimeAQ[2] << ";" << std::endl;
    }

}
