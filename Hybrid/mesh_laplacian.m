 function [lap,edge] = mesh_laplacian(vertex,face)
 
 % mesh_laplacian - Laplacian of irregular triangular mesh
 %
 % Usage: [lap,edge] = mesh_laplacian(vertex,face)
 %
 % Returns 'lap', the Laplacian (2nd spatial derivative) of an
 % irregular triangular mesh, and 'edge', the linear distances
 % between vertices of 'face'.  'lap' and 'edge' are square,
 % [Nvertices,Nvertices] in size, sparse in nature.
 %
 % It is assumed that 'vertex' contains the (x,y,z) Cartesian
 % coordinates of each vertex and that 'face' contains the
 % triangulation of vertex with indices into 'vertex' that
 % are numbered from 1:Nvertices.  For information about
 % triangulation, see 'help convhull' or 'help convhulln'.
 %
 % The neighbouring vertices of vertex 'i' is given by:
 %
 % k = find(edge(i,:));
 %
 % The math of this routine is given by:
 %
 % Oostendorp, Oosterom & Huiskamp (1989),
 % Interpolation on a triangulated 3D surface.
 % Journal of Computational Physics, 80: 331-343.
 %
 % See also, eeg_interp_scalp_mesh
 %
 
 % $Revision: 1.2 $ $Date: 2005/08/14 21:25:08 $
 
 % Licence:  GNU GPL, no implied or express warranties
 % History:  04/2002, Darren.Weber_at_radiology.ucsf.edu
 %           - initial version was inefficient and incorrect
 %             at one point.
 %           (c) 04/2002 Robert Oostenveld
 %           - completely revised/reconstructed code (lapcal.m)
 %           - agreed to release into eeg_toolbox under GNU GPL
 %           04/2002, Darren.Weber_at_radiology.ucsf.edu
 %           - modified edge initialization from sparse to
 %             full matrix and slightly improved speed of
 %             calculation for edge norms
 %           - added tic/toc timing report
 %           07/2002, Darren.Weber_at_radiology.ucsf.edu
 %           - added check for duplicate vertices (lines 78-79)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 nvertex = size(vertex,1);%the size of the first dimention, number of rows
 nface = size(face,1);%the size of the first dimention, number of rows
 
 fprintf('MESH_LAPLACIAN: Calc Laplacian matrix for %5d vertices...',nvertex);
 tic % Start stopwatch timer
 
 % the matrix 'edge' is the connectivity of all vertices
 edge = sparse(nvertex,nvertex); %generates an nvertex-by-nvertex all zero sparse matrix.
 
 for i=1:nface,
     
     % compute the length of all triangle edges (Diff is [3x3])
     Diff = [vertex(face(i,[1 2 3]),:) - vertex(face(i,[2 3 1]),:)]; 
     %face(i,[1,2,3] accesses a raw in face matrix.vertex(face(i,[1,2,3],:) 
     %accesses 3 raws in the vertex matrix. each raw is defined by 1 
     %element of the raw face(i,[1,2,3]). so this term 
     %vertex(face(i,[1,2,3]),:) returns the  Cartesian coordinates of the 
     %three vertices that form the ith face. the other term is as same as 
     %the first one but with a change in the order of raws so the 
     %difference betweeen them can give you the length of the three edges
     % that form the ith face .. the lenght in x and in y and in z.
     Norm = sqrt( sum(Diff.^2, 2) ); % square the diff matrix, sum each raw, then find the square root of the sum
     edge(face(i,1),face(i,2)) = Norm(1);% this is the first edge of the ith face
     edge(face(i,2),face(i,3)) = Norm(2);%this is the second edge of the ith face
     edge(face(i,3),face(i,1)) = Norm(3);% this is the third edge of the ith face
     % make sure that all edges are symmetric
     edge(face(i,2),face(i,1)) = Norm(1);
     edge(face(i,3),face(i,2)) = Norm(2);
     edge(face(i,1),face(i,3)) = Norm(3);
 end
 
 % Using edge to identify nearest vertices, calculate
 % the Laplacian for an irregular mesh
 lap = sparse(nvertex,nvertex);
 for i=1:nvertex,
     
     k = find(edge(i,:));        % the indices of the neighbours
     
     % remove any duplicate neighbour vertices
     [vert, m, n] = unique(vertex(k,:),'rows');
     k = k(m);
     
     ni = length(k);             % the number of neighbours
     
     hi = mean(edge(i,k));       % the average distance to the neighbours
     invhi = mean(1./edge(i,k)); % the average inverse distance to the neighbours
     
     lap(i,i) = -(4/hi) * invhi; % Laplacian of vertex itself
     
     lap(i,k) =  (4/(hi*ni)) * 1./edge(i,k); % Laplacian of direct neighbours
     
     % Laplacian is zero for all indirect neighbours
     % See Oostendorp, Oosterom & Huiskamp (1989, pp. 334-335)
 end
 
 t = toc;
fprintf('done (%6.2f sec).\n',t);
 
return
 end