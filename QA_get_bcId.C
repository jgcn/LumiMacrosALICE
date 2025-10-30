#include "CandRoot.h"
#include "GlobalVariables.h"
#include "InputFromUser.h"
#include "vdmUtilities.h"

void QA_get_bcId(int fill)
{
    // get name of files and set pointers to trees
  Set_input_file_names(fill);
  Set_pointers_to_input_files_and_trees();

  // get the ingo
  int nIBC = GetNumberInteractingBunchCrossings();
  int *bunches = new int [nIBC];
  GetBunchIndices(bunches);

  // print the info
  cout << " There are " << nIBC << " interacting bunches" << endl;
  for(int i=0;i<nIBC;i++) {
    cout << "   bcID = " << bunches[i] << endl;
  }

  // for C++ format
  cout << Form("int bunchB[%i] = {",nIBC) << endl;
  for(int i=0;i<(nIBC-1);i++) {
    cout << bunches[i] << ",";
  }
  cout << bunches[nIBC-1] << "};" <<endl;

}
