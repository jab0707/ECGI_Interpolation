clear
Training_Data=[];
badlead=130;
%% loading the geometry of the torso
load('18-08-09-192Tank.mat');
%load('tank_geom_192.mat');
vertex=tank_geom.pts';
face=tank_geom.fac';
%% loading the torso test signal and beats segments (start time and end time of each beat)
load('Run0038_5_1_1.mat')
%load('Run0017_10_1.mat');
TestSignal=Sig;
load('Run0038_5_1_1_Segments.mat');
%load('Run0017_10_1_Segments.mat');
BeatsSegments=beats;
%% generating the segmented beats
[TestBeat] = BeatSegmetation(TestSignal,BeatsSegments);
%% Generating the training set and compute the coefficient matrix (the fullset of leads) 
Run=[17 26 35 44 ];
for r=1:4
load(strcat('Run00',int2str(Run(r)),'_10_1.mat'));
Training_Data=[Training_Data Sig];
end
Training_Data(badlead,:)=Training_Data(131,:); %this is to replace the bad lead with a neighbor signal
coeff = pca(Training_Data');
%% Selecting the test electrodes
load('UtahBrokenLeads1.mat')
missed_electrodes=find(brokenLeads);
brokenLeads(badlead)=0;
Test_electrodes=find(brokenLeads);% these are the missed electrodes without the bad leads.
%% optimizing the number of eigenvectors
index=(1:192)';
[SigInt] = GenLapInterpolatedBeat(vertex,face,index,Training_Data(:,1:1000),Test_electrodes);
[eig] = OptimizingEigenNumber(SigInt,Test_electrodes,Training_Data(:,1:1000),coeff);
%% Applying Hybrid interpolation
for j=1:length(TestBeat)
[BeatLap{1,j}] = GenLapInterpolatedBeat(vertex,face,index,TestBeat{1,j},missed_electrodes);
[Sig] = HybridInterpolation(missed_electrodes,BeatLap{1,j},eig,coeff);
BeatH{1,j}=Sig;
end
%% ploting the hybrid interpolated signal vs the recorded signal
beatNum=8;
for i=1:length(missed_electrodes)
plot(TestBeat{1,beatNum}(missed_electrodes(i),:));hold on;plot(BeatH{1,beatNum}(missed_electrodes(i),:));plot(BeatLap{1,beatNum}(missed_electrodes(i),:));hold off;legend('Recorded signal','Hybrid interpolated Signal','Laplacian interpolated signal')
missed_electrodes(i)
pause
end


