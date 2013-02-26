Distance-Kernel Method (DKM)

This directory contains the Matlab code to apply distance-based clustering for model selection.   
The basic concept is to select a few representative models among a larger set via clustering (classification).  The clustering is performed based on differences in proxy responses (time serie vector) for the entire set of models. The full model (non proxy) response is evaluated on the selected medoids for uncertainty quantification.  The algorithm consists of 3 main steps:

1.	Construction of a metric space by applying multi-dimensional scaling (MDS) on the pairwise distance between models. The distance is defined as the difference in proxy responses.
2.	Kernel k-medoid clustering is then applied in the metric space and the resulting medoids are used for evaluation of the full model responses
3.	The model responses and the cluster weights are used to construct the P10-P50-P90 curves of the response(s) to be predicted. 

References: 

Scheidt C. and Caers J. (2009): Representing Spatial Uncertainty Using Distances and Kernels, Mathematical Geosciences, 41(4):397-419. DOI:10.1007/s11004-008-9186-0 

Scheidt C. and Caers J. (2009): Uncertainty Quantification in Reservoir Performance Using Distances and Kernel Methods - Application to a West-Africa Deepwater Turbidite Reservoir.  SPE J, 14 (4):680:692, (2009), SPEJ 118740-PA. DOI:10.2118/118740-PA
