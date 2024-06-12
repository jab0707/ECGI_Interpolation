function [SigH] = HybridInterpolation(missed_electrodes,SigLap,eig,coeff)
dataAfterLap=SigLap;
PCs  = dataAfterLap' *coeff;
reconMap = (PCs(:,1:eig)) * coeff(:,1:eig)';
reconMap=reconMap';
dataAfterLap(missed_electrodes,:)=reconMap(missed_electrodes,:);
SigH=dataAfterLap;
end

