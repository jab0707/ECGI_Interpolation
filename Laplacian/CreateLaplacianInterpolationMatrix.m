function [lapInterp] = CreateLaplacianInterpolationMatrix(geom,badLeads,pathLength,lap)
%CREATELAPLACIANINTERPOLATIONMATRIX Summary of this function goes here
%   Detailed explanation goes here

if isfield(geom,'pts')&&~isfield(geom,'node')
    geom.node = geom.pts;
end
if isfield(geom,'fac')&&~isfield(geom,'face')
    geom.face = geom.fac;
end

if size(geom.node,1) ~= 3
    if size(geom.node,2) == 3
        geom.node = geom.node';
    else
        throw('No dimensions of the input geometry correspond to [x y z]\n')
    end
end
if size(geom.face,1) ~= 3
    if size(geom.face,2) == 3
        geom.face = geom.face';
    else
        throw('No dimensions of the input geometry faces are as expected\n')
    end
end
if ~exist('lap')
    adj = computeAdjacencyMatrix(geom,pathLength);

    weightFunct = @(index) adj(index,:);
    [D,H] = meshVolDiffHessMatrix(geom,weightFunct);
    lapMatrix = LaplacianMatrixFromHessianMatrix(H);

    
else
    lapMatrix = lap;
end
known = find(badLeads == 0);
unknown = find(badLeads == 1);
L_k = lapMatrix(known,known);
L_u = lapMatrix(unknown,unknown);
L_c = lapMatrix(known,unknown);

%maps from known to unknown
interpMat = (L_c'*L_c + L_u'*L_u)\(-(L_c'*L_k + L_u'*L_c'));

lapInterp = zeros(size(geom.node,2),1);
lapInterp(known) = 1;
lapInterp = diag(lapInterp);
lapInterp(unknown,known) = interpMat;
lapInterp(:,unknown) = [];
end

