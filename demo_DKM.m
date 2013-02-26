%% Distance-Based Methods for Uncertainty Quantification
% Author: Celine Scheidt
% Date: August 2011
% Updated: July 2012


% Below is where you load the data, it must be in this format:
% - ProxyResponse:  Matrix containing the proxy responses for each model. One row is
%                     for one response. The size is  Number of models (n) x Number of timesteps for proxy
% - FullResponse: Matrix containing the full responses for each model. One row is
%                     one response. The size is  Number of models (n) x
%                     Number of timesteps for full response
% - TimeF: Timesteps for the full responses. Used only to display the quantiles in the provided example.


load 'InputData_hydro.mat'  % Hydro example with tracer concentrations
%load 'InputData_WCA.mat'   % Oil example, with field cumulative oil production


 % Compute the distance using Proxy Responses
d = pdist(ProxyResponse);

% Multi-Dimensional Scaling (MDS) using distance d
[Xd_, e_] = cmdscale(d);

% Reduction of the dimension of the MDS space, here 3D
dims = 3;
Xd = Xd_(:,1:dims);

% Definition of the kernel matrix - rbf Gaussian
sigma = std(reshape(d,[],1));  % bandwidth
K = rbf_kernel(Xd,sigma);

% Perform the clustering (kernel k-medoid)
nbclusters = 10;  % number of clusters to create (7 for WCA)
Clustering = kkmedoid(Xd,nbclusters,25,K);

% Plot of the clusters and medoids
plotcmdmap(Xd, Clustering)


% Computation of the quantiles (P10-P50-P90) for the full set (as a reference) and for the selected
% responses
Quantiles_ref = QuantileComputation(FullResponse,[0.1,0.5,0.9]); 
Quantiles_est = QuantileComputation(FullResponse,[0.1,0.5,0.9],Clustering);

figure; axes('FontSize',12);hold on;
h1 = plot(TimeF,Quantiles_ref,'-.k','LineWidth',4);
h2 = plot(TimeF,Quantiles_est ,'-r','LineWidth',3);
legend([h1(1),h2(1)],'Initial Set','Selected Responses','location','NorthWest')

