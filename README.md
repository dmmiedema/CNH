# CNH
Method to infer Copy Number intra-tumor Heterogeneity from a single chromosomal copy number measurement.

# LICENSE
Please read the end user license agreement before usage.

# DESCRIPTION
MATLAB function to calculate copy number intra-tumor heterogeneity (CNH) from a single copy number measurement.
In this method, relative segmented copy numbers are transformed to absolute copy numbers and the distance
to integer values is measured for each segment. This procedure is
performed for either a range of ploidies/purities or for fixed
ploidy/purity, in case either or both quantities are known from previous
measurements.
Copy number heterogeneity (CNH) is defined as the minimum average distance of
segments to the closest integer value, weighted by segment length. 
In addition to CNH, this function returns the ploidy and purity
corresponding to the inferred CNH.

# INPUT
1st input argument:   (seg_val)       values of relative copy numbers per segment, vector of size Nx1
2nd input argument:   (seg_len)       segment lengths, vector of size Nx1
3th input argument:   (ploidy)        tumour ploidy, use empty vector ( = [] ) for grid search 
4th input argument:   (purity)        sample purity, use empty vector ( = [] ) for grid search     

# OUTPUT
1st output argument:  (CNH_out)       inferred CNH 
2nd output argument:  (ploidy_out)    inferred ploidy for empty input ploidy, otherwise same as input ploidy. 
3th output argument:  (purity_out)    inferred purity for empty input purity, otherwise same as input purity.

