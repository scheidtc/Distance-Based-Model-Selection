%
% This function applies a kernel k-medoid algorithm (clustering) 
% in the Feature Space defined by the kernel matrix K. 
 
% Author: Celine Scheidt
% Date: April 2009
% Updated: July 2012


function Clustering = kernel_kmedoid(Xd,nbclusters,K)

%% Input Parameters
%   - Xd: location of the points to be plotted (from MDS). Rows of Xd are 
%         the coordinates in p-dimensional space. Dimension: Number of
%         realizations (n) x p
%   - nbclusters: number of clusters to construct
%   - K: Kernel matrix 

%% Output Parameters 
%   - Clustering: Results of the clustering, which contains:
%        - T: is a vector of length n, which contains the cluster index that each model belongs to
%        - medoids: is a vector of length nbcluster containing the index
%                   of the medoids
%        - weights: is a vector of ength nbcluster containing the number
%                   of models in each cluster.
%


maxIterations = 50;
npoints = size(Xd,1);
totsumDBest = Inf;

for iter = 1:100 % 100 clustering are performed, the best is returned

    
    %% 1. Random selection of the inital medoids.

    initMedoids = randperm(npoints);
    initMedoids = initMedoids(1:nbclusters);


    %% 2. Associate each points to the closest medoid

    % 2.1 Compute distance between each point and the selected initial medoids. Note that the distance is in the 
    % feature space, so it can be computed using the kernel matrix only 

    dist_points_medoids = dist_feature(K,1:npoints,initMedoids);

    % 2.2 Minimization: T contains the cluster index each model belongs to
    [~,T] = min(dist_points_medoids,[],2);


    %% 3.  Update the medoid of each cluster and re-assign each point to a cluster

    Medoids_prev_iter = ones(1,nbclusters);
    currentMedoids = initMedoids;
    nbIter = 0;

    while (~all(currentMedoids == Medoids_prev_iter) && nbIter < maxIterations)  % while cluster configuration is changing and maxIteration not reached
        
        Medoids_prev_iter = currentMedoids;

        % For each cluster
        for i = 1:nbclusters
            pts_in_cluster = find(T == i);
            
            % Compute distance between each point in the cluster and its medoid.
            dist_within_cluster = dist_feature(K,pts_in_cluster,pts_in_cluster);
            
            % minimize the distance and select the new medoid
            [~,idx_min] = min(mean(dist_within_cluster));
            currentMedoids(i) = pts_in_cluster(idx_min);
        end

        % New medoids are defined, compute distances between the points and new medoids
        dist_points_medoids = dist_feature(K,1:npoints,initMedoids);
        
        %Associate each point to the closest medoid
        [~,T] = min(dist_points_medoids,[],2);
        nbIter = nbIter +1;    

    end

    % Compute the sum of distances
    sumD = zeros(nbclusters,1);
    for i = 1:nbclusters
        sumD(i) = sum(dist_points_medoids(T == i,i));
    end
    totsumD = sum(sumD);

    if totsumD < totsumDBest  % Keep the best clustering so far
        totsumDBest = totsumD;
        TBest = T;
        currentMedoidsBest = currentMedoids;
    end

end
%% Once the medoids are defined, store the outputs

weights = zeros(nbclusters,1);
for i = 1:nbclusters
    weights(i) = sum(TBest == i);
end
    

Clustering.T = TBest;  
Clustering.medoids = currentMedoidsBest;  
Clustering.weights = weights;

end

%% This function compute distance in Feature Space using the Kernel function
% d_f(x,y) = K(x,x) -2K(x,y) + K(y,y)

function d_f = dist_feature(K,x,y)  % x and y can be vectors

    k1 = diag(K(x,x));
    k2 = diag(K(y,y));
    d_f = k1(:,ones(1,length(y))) -2*K(x,y) + k2(:,ones(length(x),1))';

end

