function [eig] = OptimizingEigenNumber(BeatAfterLap,Test_electrodes,TrainingBeat,coeff)
% TrainingBeat: the first beat that we used to calculate the coeff
% dataAfterLap: the laplacian interpolated signals of the first beat
% Test_electrodes: electrodes that we want to measure the performance of
% our method on them to  optimize the number of eigenvectors.

PCs  = BeatAfterLap' *coeff;
for eig = 1:20
    reconMap = (PCs(:,1:eig)) * coeff(:,1:eig)';
    reconMap=reconMap';
    BeatAfterLap(Test_electrodes,:)=reconMap(Test_electrodes,:);
    LowerNumSamples=min(size(TrainingBeat,2),size(BeatAfterLap,2));
    for m=1:length(Test_electrodes)
        R=corrcoef(TrainingBeat(Test_electrodes(m),1:LowerNumSamples),BeatAfterLap(Test_electrodes(m),1:LowerNumSamples));
        R=R(1,2);
        CorrOnebeat(m)=R;
    end
    CC(eig)=median(CorrOnebeat);
end
CC(1:4)=0;[M1,I1]=max(CC);
eig=I1;
end



