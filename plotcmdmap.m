function plotcmdmap(Xd,KKM)

%% Function which maps the points defined in MDS (Xd) and (optional) colors the
%% points according to the clustering results obtained by kernel_kmedoid.m

% Author: Celine Scheidt
% Date: April 2009


%% Input Parameters:
% - Xd: location of the points to be plotted (from MDS). Rows of Xd are the coordinates in p-dimensional space
% - KKM: results from the clustering in the format given by the function kernel_kmedoid. This is optional.  
    % If given, point are colored by clusters; medoids are represented by squares 

if nargin ==2
    
    % definition of the color of each cluster
    C  = colormap(jet(64));
    C = C(floor(linspace(1,64,max(KKM.T))),:);

    figure
    for i=1:length(KKM.medoids)
        plot(Xd(KKM.T==i,1), Xd(KKM.T==i,2), 'o', 'Color',C(i,:), ...  % points in clusters
            'MarkerSize', 8, 'LineWidth', 1.3,'MarkerFaceColor',C(i,:), ...
            'MarkerEdgeColor','k');

        hold on
        plot(Xd(KKM.medoids(i),1), Xd(KKM.medoids(i),2), 'bs', ...          % medoid
            'MarkerSize', 12, 'LineWidth', 2,'MarkerFaceColor',C(i,:), ...
            'MarkerEdgeColor','k');
    end

else
        
    figure
    plot(Xd(:,1), Xd(:,2), 'o', 'MarkerSize', 8, 'LineWidth', 1.2, ...
        'MarkerFaceColor',[0.7 0.7 0.7], 'MarkerEdgeColor','k');

end


hold off
grid off

set(gca,'LineWidth',3)
set(gca,'XTickLabel',{''})
set(gca,'YTickLabel',{''})

end
