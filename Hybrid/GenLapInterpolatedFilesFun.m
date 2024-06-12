function [SigInt, All_missed_electrodes] = GenLapInterpolatedFilesFun(vertex,face,FiltSig,Test_electrodes,missed_electrodes)
[lap,edge] = mesh_laplacian(vertex,face);
index=(1:192)';
index_after_del=index;
data_aft_del=FiltSig;
All_missed_electrodes=[missed_electrodes;Test_electrodes];
index_after_del(All_missed_electrodes)=[];
[int, keepindex, repindex] = mesh_laplacian_interp(lap, index_after_del);
data_aft_del(All_missed_electrodes,:)=[];
SigInt=int*data_aft_del;
end

