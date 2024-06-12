function [SigInt] = GenLapInterpolatedBeat(vertex,face,index,beat,missed_electrodes)
[lap,edge] = mesh_laplacian(vertex,face);
index_after_del=index;
data_aft_del=beat;
index_after_del(missed_electrodes)=[];
[int, keepindex, repindex] = mesh_laplacian_interp(lap, index_after_del);
data_aft_del(missed_electrodes,:)=[];
SigInt=int*data_aft_del;
end

