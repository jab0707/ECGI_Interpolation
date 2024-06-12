function [] = GeneratingHybridIntFiles(OriginalFileName,Run,Beat,savefiles)
%% loading the geometry of the torso
load('D:\Collaboration for Interpolating low amplitude signal\Utah1 Data\tank_geom_192.mat');
vertex=tank_geom.pts;
face=tank_geom.fac;
BeatH=cell(1,length(Beat));
coeff=cell(1,length(Beat));
%% compute the coefficient matrix using all the signals as training set 1st method
Training_Data=Beat{1,1};
Training_Data(130,:)=Training_Data(131,:);
coefficient = pca(Training_Data');
for i=1:length(Beat)
coeff{1,i}=coefficient;
end

%% Initiate a loop for BrokenType   
for L=50%[1:5:190]
%% Selecting the test electrodes
% This will help in optimizing the number of eigenvectors that should be selected
load(strcat('D:\Collaboration for Interpolating low amplitude signal\Utah2 Data\Run00',num2str(Run),'\Low signals indices files\LowSig',int2str(L),'.mat'))
missed_electrodes=LowSig';
missed_electrodes_No130=missed_electrodes(1:end-1);
[Test_electrodes] = SelectingTestElectrodes(missed_electrodes_No130,face);
%% optimizing the number of eigenvectors
Devisions=floor(length(Test_electrodes)/5);
Test_electrodes_Devisions = reshape(Test_electrodes(1:(5*Devisions)),[Devisions,5]);
for j=1:length(Beat)
for i=1:5
 [SigInt, All_missed_electrodes] = GenLapInterpolatedFilesFun(vertex,face,Beat{1,j},Test_electrodes,missed_electrodes);
 [eig(i),RMSE_All] = EigNumOptimizing(SigInt, All_missed_electrodes,Test_electrodes_Devisions(:,i),Beat{1,j},coeff{1,j});
end
eigAve=round(sum(eig)/5);
if eigAve==0
    eigAve=1;
end
%% Applying Hybrid interpolation
[SigLap, All_missed_electrodes2] = GenLapInterpolatedFilesFun(vertex,face,Beat{1,j},[],missed_electrodes);
[Sig] = HybridInterpolation(missed_electrodes,SigLap,eigAve,coeff{1,j});
BeatH{1,j}=Sig;
end

if savefiles==1
save(strcat('D:\Collaboration for Interpolating low amplitude signal\Utah2 Data\Run00',int2str(Run),'\1.Torso Signals\Hybrid interpolated files\',OriginalFileName(1:end-4),'_H',int2str(L)),'BeatH');    
end
    
end

end

