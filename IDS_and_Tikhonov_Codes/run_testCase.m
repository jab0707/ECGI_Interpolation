% Sample script to run IDS interpolation and Tikhonov regularization
%
% Yesim Serinagaoglu (yserin@metu.edu.tr)
% 25/06/2024

close all
clear all

% Load sample data
load("sampleData.mat") % interpInputs
%       geomTorso: [1×1 struct]         torso geometry (nodes and triangles) mm
%              Fs: 2048                 sampling freq
%       exclLeads: 44                   global bad lead (already included in the BL)
%             Dep: [40 450]             QRS range
%             Rep: [450 840]            T-wave range
%        BaseLine: [1000 1021]          baseline range
%         bspBeat: [128×1128 double]    single beat 128-lead BSP
%     brokenLeads: [128×1 double]       BL3 scenario broken leads
%         forwMat: [128×894 double]     forward matrix


% Generate reduced BSP for the good leads
badLeadIndices = find(interpInputs.brokenLeads==1);
goodLeadIndices = find(interpInputs.brokenLeads==0);
reducedData = interpInputs.bspBeat(goodLeadIndices,:);

% Apply IDS interpolation
choi = 'Near'; % OR use 'All'
switch choi
    case 'Near'
        Radius = 10; % mm
        interpolatedBSP = interpolationIDS(reducedData,interpInputs.geomTorso.node,interpInputs.brokenLeads, Radius);
    case 'All'
        interpolatedBSP = interpolationIDS(reducedData,interpInputs.geomTorso.node,interpInputs.brokenLeads);
end

% Plots
for k=1:length(badLeadIndices)
    plot(interpInputs.bspBeat(badLeadIndices(k),:),'k');
    hold on
    plot(interpolatedBSP(badLeadIndices(k),:),'r');
    hold off
    pause
end

%%
% Inverse sln

% data prep
choi = 'IDS'; % Alternatives: 'RMV' and 'IDS'
switch choi
    case 'FULL'
        % interpInputs.exclLeads is measured bad lead and should be excluded
        A = interpInputs.forwMat;
        A(interpInputs.exclLeads,:) = [];
        Sig = interpInputs.bspBeat;
        Sig(interpInputs.exclLeads,:) = [];
    case 'RMV'
        % remove all broken leads - reduce forward and data matrices
        A = interpInputs.forwMat(goodLeadIndices,:);
        Sig = interpInputs.bspBeat(goodLeadIndices,:);
    case 'IDS'
        A = interpInputs.forwMat;
        Sig = interpolatedBSP;
end

% Run Tikhonov regularization
addpath('RegTools_v4\')
[X_varLambda, X_medLambda, allLambdaVals, medianLambdaVal] = tikhonovRegTools(Sig, A); 
% You may also specify a frame range for taking the median value of all lambdas as an input
