clear
%% loading the files of tank1 geometry and data
% load('case2_geom_128.mat')%for case2
% load('tank2_geom_192.mat')% for tank2 
% load('Bordeaux_data_segments.mat')%for case2
% load('Tank2_beat_segmentation_new.mat')% for tank2
load('Bordeaux_tank3_beat_segments_corrected.mat')%for tank4
load('case4_geom_128.mat')%for tank4
% badlead=130;% for tank2 
badlead=44;%for case2 & case 4
% FileName='Run0038-ts_';% for tank2 
FileName='Case4_Pacing_';%for case4
% BrokenFileName='UtahBrokenLeads';% for tank2 
BrokenFileName='BordeauxBrokenLeads';%for case4
% TankNum='tank2';% for tank2 
TankNum='case4';%for case2
vertex=tank_geom.pts;
face=tank_geom.fac;
% Intiate a loop for filter type
for FiltType=5:5
% compute the coefficient matrix using all the signals as training set 1st method
load(strcat(FileName,num2str(FiltType),'_1_1'));
Training_Data=Sig(:,beats{1, 1}.fromIdx:beats{1, 1}.toIdx);
% Training_Data(badlead,:)=Training_Data(131,:);%for tank2
Training_Data(badlead,:)=Training_Data(111,:);%for case4
coeff = pca(Training_Data');
% Initiate a loop for BrokenType   
for BrokenType=2:7
%% Selecting the test electrodes
% This will help in optimizing the number of eigenvectors that should be selected
load(strcat(BrokenFileName,int2str(BrokenType-1)),'brokenLeads')
missed_electrodes=find(brokenLeads);
brokenLeads(badlead)=0;
Test_electrodes=find(brokenLeads);
load(strcat(FileName,num2str(FiltType),'_',int2str(BrokenType),'_2.mat'))
% this is to generate the part of the signal on the baseline
SigLap=ts.potvals;
[Sig] = HybridInterpolation(missed_electrodes,SigLap,7,coeff);
for j=1:length(beats)
SigInt=SigLap(:,beats{1, j}.fromIdx:beats{1, j}.toIdx);
OriginalSignal=Training_Data;
[eig] = OptimizingEigNum(SigInt,Test_electrodes,OriginalSignal,coeff);
%% Applying Hybrid interpolation
[BeatH] = HybridInterpolation(missed_electrodes,SigInt,eig,coeff);
Sig(:,beats{1, j}.fromIdx:beats{1, j}.toIdx)=BeatH;
end
save(strcat('C:\Users\Ali\Desktop\interpolation\Utah Data for Interpolation Group\Hybrid Interpolated Files\',TankNum,'\',FileName,int2str(FiltType),'_',int2str(BrokenType),'_6'),'Sig')
end
end