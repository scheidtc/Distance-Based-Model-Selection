%
% This function computes the quantiles of Response as a function of time.
% If clustering is performed, each response is weighted by the number of
% points in the corresponding cluster.

% Author: Celine Scheidt
% Date: April 2009


function Quantiles = QuantileComputation(Response,p,Clustering)

%% Input Parameters
%   - Response: matrix of the responses as a function of time. One line is one realization.
%   - p: scalar or vector of cumulative probability values
%   - Clustering: results of the clustering (as given by function kernel_kmedoid). Optional.
%                 If not given, return the quantiles of Response without ponderations by the clusters weights.  

%% Output Parameters 
%   - Quantiles: quantiles of the values in response.


if nargin == 2 % No clustering: each model as a weight of one
    Quantiles=quantile(Response,p);
    
else % Each model is weighted by the number of point in the cluster (weights)
    WeigthedResponse = [];
    for i = 1:size(Clustering.medoids,2)
        WeigthedResponse=vertcat(WeigthedResponse,repmat(Response(Clustering.medoids(i),:),Clustering.weights(i),1));
    end
    
    Quantiles=quantile(WeigthedResponse,p);    
end
    
end