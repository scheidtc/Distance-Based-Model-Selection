%% Distance-Based Methods for Uncertainty Quantification
% Author: Celine Scheidt
% Date: August 2011
% Updated: July 2012


% Below is where you load the data, it must be in this format:
% - SimulationProxy:  Matrix containing the proxy simulations. One row is
%                     one realization. The size is  Number of realizations x Number of timesteps for proxy
% - SimulationFull: Matrix containing the full simulations. One row is
%                     one realization. The size is  Number of realizations x Number of timesteps
%                     for fine simulations. 
% - TimeF: Timesteps for the fine simulation (to plot the quantiles)


load 'InputData.mat'

 % Compute the distance using Proxy Simulations
d = pdist(SimulationProxy);

% Multi-Dimensional Scaling (MDS) using distance d
[Xd_, e_] = cmdscale(d);

% Reduction of the dimension of the MDS space, here 3D
dims = 3;
Xd = Xd_(:,1:dims);

% Definition of the kernel matrix - rbf Gaussian
sigma = std(reshape(d,[],1));  % bandwidth
K = rbf_kernel(Xd,sigma);

% Perform the clustering (kernel k-medoid)
nbclusters = 10;  % number of clusters to create
Clustering = kernel_kmedoid(Xd,nbclusters,K);

% Plot of the clusters and medoids
plotcmdmap(Xd, Clustering)


% Computation of the Quantiles for the full set (as a reference) and for the selected
% realizations
Quantiles_ref = QuantileComputation(SimulationFull,[0.1,0.5,0.9]);
Quantiles_est = QuantileComputation(SimulationFull,[0.1,0.5,0.9],Clustering);

figure; hold on
h1 = plot(TimeF,Quantiles_ref,'-.k','LineWidth',4);
h2 = plot(TimeF,Quantiles_est ,'-r','LineWidth',3);
legend([h1(1),h2(1)],'Initial Set','Selected Responses')

