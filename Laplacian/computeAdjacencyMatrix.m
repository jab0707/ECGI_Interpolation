%% HELP:
%
%	This functionl computes the adjacency matrix on the graph defined in
%	the geometry "geom".
%	The adjacency matrix can be extended to all the neighbors of a node
%	that are at a distance "pathLength" from that node.
%
%
%	INPUT:
%		- geom - struct - input geometry. Assumed to be a triangular
%		surface mesh. The struc is composed of two fields: node and face.
%		- pathLength - int - number of jumps on the graph that define the
%		adjacency.
%
%	OUTPUT:
%		- AdjMtrx - <M,M>bool - matrix defining which nodes are neighbors
%		of which.
%
%	AUTHOR:
%		Jaume Coll-Font <jcollfont@gmail.com>
%       Updated: Jake Bergquist <jake.a.bergquist@gmail.com>
%

function [adjMtrx] = computeAdjacencyMatrix(geom, pathLength)
    
	%% define
        if ~exist('pathLength','var')
            pathLength = 1;
        end
        [geom,nodenames,elemNames] = checkMeshNodeFaceBasis_jb(geom,'','row');
        geom.node = geom.(nodenames{1});
        if ~strcmpi(elemNames{1},'face')
            if size(geom.(elemNames{1}),1)~=3
                error('Faces not found')
            end
            geom.face = geom.(elemNames{1});
        end
		[d,M] = size(geom.node);
		[~,numFac] = size(geom.face);
			
		if d ~= 3
			numFac = d;
			geom.face = geom.face';
			geom.node = geom.node';
		end
			
	%% create 1rst order adjacency matrix
		adjMtrx = false(M);
        adjMtrx(sub2ind([M,M],geom.face(1,:),geom.face(2,:))) = 1;
        adjMtrx(sub2ind([M,M],geom.face(1,:),geom.face(3,:))) = 1;
        adjMtrx(sub2ind([M,M],geom.face(2,:),geom.face(1,:))) = 1;
        adjMtrx(sub2ind([M,M],geom.face(2,:),geom.face(3,:))) = 1;
        adjMtrx(sub2ind([M,M],geom.face(3,:),geom.face(1,:))) = 1;
        adjMtrx(sub2ind([M,M],geom.face(3,:),geom.face(2,:))) = 1;
% 		for ii = 1:numFac;
% 			adjMtrx(geom.face(1,ii),geom.face(2,ii)) = 1;%distance(geom.face(1,ii),geom.face(2,ii));
% 			adjMtrx(geom.face(2,ii),geom.face(3,ii)) = 1;%distance(geom.face(2,ii),geom.face(3,ii));
% 			adjMtrx(geom.face(3,ii),geom.face(1,ii)) = 1;%distance(geom.face(3,ii),geom.face(1,ii));
% 			adjMtrx(geom.face(2,ii),geom.face(1,ii)) = 1;%distance(geom.face(2,ii),geom.face(1,ii));
% 			adjMtrx(geom.face(3,ii),geom.face(2,ii)) = 1;%distance(geom.face(3,ii),geom.face(2,ii));
% 			adjMtrx(geom.face(1,ii),geom.face(3,ii)) = 1;%distance(geom.face(1,ii),geom.face(3,ii));
% 		end
		
	%% extend adjacency matrix according to path length
	%AdjMtrx = adjMtrx + eye(M);
    %AdjMtrx = adjMtrx;
    adjMtrx(1:M+1:M*M) = 1;
	for rep = 1:(pathLength-1)
		for n =1:M
			neighbors = find(adjMtrx(n,:));
			for m = neighbors
				adjMtrx(n,:) = adjMtrx(n,:)|adjMtrx(m,:);
			end
		end
	end
end