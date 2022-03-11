# ZPrimeBoosted

Codes are in folder 'codes'. Folder 'run' contains all the outputs

  1) the dataset is read and the selection is applied in 'run_Selection.py'. 
     Here each tree contained within the corresponding root file is converted into a RDataFrame, then all the selections are applied sequentially and finally the        visible di-top invariant mass is constructed and the histogram is saved in a root file.
     
 2)  the Extended ML fit for several mass points signal strenght is implemented in 'MaximumLikelihoodFit.py'. Here workspaces are created and saved in a root file.

 3)  The Hypothesis test for Discovery (Background-only hypothesis) and for upper limit setting on signal strenght (Signal + Background hypothesis, for several signal mass values) is performed with the 'StatisticalTests.py' code.
