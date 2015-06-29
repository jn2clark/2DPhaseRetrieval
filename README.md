# 2DPhaseRetrieval Jesse Clark
2D Phase Retrieval code to reconstruct images from x-ray diffraction data.

Uses iterative phase retrieval based on projections onto constraint sets (the 'support' and 'modulus' constraint).  This version is suitable for 2D X-ray free electron laser data with a missing central stop.

Also implements some global optimization via guided methods.  See MatlabPhasing_v2_1.m for specific details of the parameters.
  
To run, set the Matlab path and then run MatlabPhasing_v2_1.m with default parameters.
