function [ normals, iVect, jVect ] = normalVectorsJD( tri, method, neighs )
%NORMALVECTORSJD calculate vertex normals on the mesh
% Inputs: 
%       tri    : TriRep of mesh
%       method : Method of vertex normal approximation (optional)
% 
% Modified by JD 18/11/14
    
    % Default method is non-weighted norm averaging
    if nargin < 2 || isempty(method) ; method = ''; end
    if (nargin < 3 || isempty(neighs)) && nargout > 1; neighs = neighborsJD(tri); end
    nFace = size(tri.Triangulation, 1);
    nVert = size(tri.X, 1);
    normals = zeros(nVert,3);

    for i=1:nFace
        f = tri.Triangulation(i,:);
        % Compute the normal to the face
        edge1 = tri.X(f(3),:)-tri.X(f(1),:);
        edge2 = tri.X(f(2),:)-tri.X(f(1),:);
        n = cross(edge1, edge2);
        if strcmp(method, 'weighted')
            squares = (edge1.^2) * (edge2.^2)';
            n = n/squares;
        else
            n = n/norm(n);
        end
        
        % Sum calculated vector for all adjacent vertices
        for j=1:3
            normals( f(j),: ) = normals( f(j),: ) + n;
        end
    end

    % Normalize
    for i=1:nVert
    n = normals(i,:);
    normals(i,:) = n/norm(n);
    end
    
    % Check vectors point outwards (via sum of dot product vith rays)
    rays = tri.X - repmat(mean(tri.X, 1), nVert, 1);
    sumdot = sum(dot(normals, rays, 2), 1);
    if sumdot < 0, normals = -normals; end
    
    % Compute orientation vectors for each point
    if nargout > 1
        iVect = nan(nVert, 3);  jVect = nan(nVert, 3);
        hasneighs = ~cellfun(@isempty, neighs)';
        for i = find(hasneighs)
            tmpVect = tri.X(neighs{i}(1),:) - tri.X(i,:);
            jVect(i,:) = cross(tmpVect, normals(i,:));
            jVect(i,:) = jVect(i,:) ./ norm(jVect(i,:));
            iVect(i,:) = cross(normals(i,:), jVect(i,:));
        end
    end
end

