function [transferA, virtualA, virtualG, virtual3DG, virtualL, virtual3DL] = compute(obj)
%COMPUTE method - calculate transfer matrix from input data
%formulations : standard or method of fundamental solutions
%

%%  METHOD OF FUNDAMENTAL SOLUTIONS (Wang & Rudy '06)
torsoX = obj.pciCase.triTorso.X(obj.pchannelMask,:);
torsoNormals = obj.pciCase.normTorso(obj.pchannelMask,:);
heartX = obj.pciCase.triHeart.X;
heartNormals = obj.pciCase.normHeart;
if all(obj.pchannelMask)
    torsoTri = obj.pciCase.triTorso.Triangulation;
else
    torsoTri = convexHull(delaunayTriangulation(torsoX));
end

nptHeart = size(heartX, 1); nptTorso = size(torsoX, 1);

% Calculate center of heart mesh for shrunk/inflated meshes
meanX = mean(heartX, 1);
meanTX = mean(torsoX, 1);

% Inflate torso mesh and contract heart mesh for ficticious sources
nptSources = nptHeart + nptTorso;
heartRatio = 0.8;   torsoRatio = 1.2;
fictHeart = repmat(meanX, nptHeart, 1) + ...
    heartRatio * (heartX - repmat(meanX, nptHeart, 1));
fictTorso = repmat(meanX, nptTorso, 1) + ...
    torsoRatio * (torsoX - repmat(meanX, nptTorso, 1));
fictSources = [fictHeart ; fictTorso];

% Calculate matrix A (equation Aa = b)
% 1 ) Equations corresponding to Dirlichet condition
drlchSqTorso = repmat(sum(torsoX .^ 2, 2),1,nptSources);
drlchSqSource = repmat(sum(fictSources' .^ 2, 1),nptTorso,1);
drlchSourceTorso = torsoX * fictSources';
drlchA = 1./sqrt(drlchSqTorso + drlchSqSource - 2 * drlchSourceTorso);

% 2 ) Equations corresponding to Neumann condition
nmnSub1 = torsoNormals * fictSources';
nmnSub2 = repmat(sum(torsoX .* torsoNormals, 2),1,nptSources);
nmnA = (nmnSub1-nmnSub2).*drlchA.^3;

% 3 ) Pile up and add vector for a0
transferA = [[ ones(nptTorso, 1) drlchA ];...
    [zeros(nptTorso, 1)  nmnA  ]];

% Compute matrix from virtual sources to epicardium
forSqHeart = repmat(sum(heartX.^2,2),1,nptSources);
forSqSource = repmat(sum(fictSources'.^2,1),nptHeart,1);
forSourceHeart = heartX*fictSources';
forA = 1./sqrt(forSqHeart+forSqSource - 2*forSourceHeart);
virtualA = [ones(nptHeart, 1) forA];
assignin('base','virtualA',virtualA);
assignin('base', 'heartX', heartX);
assignin('base', 'fictSources', fictSources);


forDx = repmat([1 0 0], nptHeart, 1) * fictSources';
forDx2 = repmat(sum(heartX.* repmat([1 0 0], nptHeart, 1), 2),1,nptSources);

forDy = repmat([0 1 0], nptHeart, 1) * fictSources';
forDy2 = repmat(sum(heartX.* repmat([0 1 0], nptHeart, 1), 2),1,nptSources);

forDz = repmat([0 0 1], nptHeart, 1) * fictSources';
forDz2 = repmat(sum(heartX.* repmat([0 0 1], nptHeart, 1), 2),1,nptSources);

forDNorm = heartNormals * fictSources';
forDNorm2 = repmat(sum(heartX.* heartNormals, 2),1,nptSources);

% Compute matrix from virtual sources to gradient/laplacian
virtual3DG(:,:,1) = [zeros(nptHeart, 1) -(forDx-forDx2).*forA.^3];
virtual3DG(:,:,2) = [zeros(nptHeart, 1) -(forDy-forDy2).*forA.^3];
virtual3DG(:,:,3) = [zeros(nptHeart, 1) -(forDz-forDz2).*forA.^3];
virtualG = [zeros(nptHeart, 1) -(forDNorm-forDNorm2).*forA.^3];

virtual3DL(:,:,1) = [zeros(nptHeart, 1) 3*((forDx-forDx2).^2).*forA.^5];
virtual3DL(:,:,2) = [zeros(nptHeart, 1) 3*((forDy-forDy2).^2).*forA.^5];
virtual3DL(:,:,3) = [zeros(nptHeart, 1) 3*((forDz-forDz2).^2).*forA.^5];
virtualL = [zeros(nptHeart, 1) 3*((forDNorm-forDNorm2).^2).*forA.^5];

%transferA = NaN;


end


