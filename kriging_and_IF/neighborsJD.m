function [ neighbors, distneighs ] = neighborsJD(tri, edgemask)
%NEIGHBORSJD Calculate neighbors for all points on the mesh
%   In case of a remeshed surface (artificially closed), neighbors will be
%   removed if they correspond to an artificial closure (ie. valves) as
%   specified by the edgemask input variable
%
% Inputs:
%       tri: input mesh
%       edgemask: input edge mask (for remeshing/artificial deconnection between points)
%
% JD 09/12/14

    if nargin < 2, edgemask = []; end
    edgemask = sort(edgemask, 2);
    
    neighbors = cell(size(tri.X, 1), 1);
    if nargout > 1, distneighs = cell(size(tri.X, 1), 1); end
    for i = 1:size(tri.Triangulation, 1)
        vertices = sort(tri.Triangulation(i, :), 2);
        
        % List pairs of points
        points = [[vertices(1) vertices(2)]; ...
                  [vertices(2) vertices(3)]; ...
                  [vertices(1) vertices(3)]];
        
        % Remove those from the mask if necessary
        if ~isempty(edgemask), points = setdiff(points, edgemask, 'rows'); end
        
        % Add the ones left to the neighbors list
        for j = 1:size(points, 1)
            neighbors{points(j, 1)} = [neighbors{points(j, 1)} points(j, 2)];
            neighbors{points(j, 2)} = [neighbors{points(j, 2)} points(j, 1)];
        end
        
        % Compute distances if necessary
        if nargout > 1
            for j = 1:size(points, 1)
                dist = tri.X(points(j, 1), :) - tri.X(points(j, 2), :);
                dist = sqrt(sum(dist.^2));
                distneighs{points(j, 1)} = [distneighs{points(j, 1)} dist];
                distneighs{points(j, 2)} = [distneighs{points(j, 2)} dist];
            end
        end
        
    end
    if nargout <= 1
        neighbors = cellfun(@unique, neighbors, 'UniformOutput', 0);
    else
        [neighbors, uia] = cellfun(@unique, neighbors, 'UniformOutput', 0);
        distneighs = cellfun(@(x,y) x(y), distneighs, uia, 'UniformOutput', 0);
    end    
end
