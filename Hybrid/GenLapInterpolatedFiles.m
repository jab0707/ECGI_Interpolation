clear
load('tank_geom_192.mat');
vertex=tank_geom.pts;
face=tank_geom.fac;
[lap,edge] = mesh_laplacian(vertex,face);
index=(1:192)';
% initiate a loop for type of filtering
for aa=5
    load(strcat('Run0053-ts_',num2str(aa),'_1_1'));
%initiate a loop for type of broken leads    
for ii=2:7
    index_after_del=index;
    data_aft_del=Sig;
    load(strcat('UtahBrokenLeads',num2str(ii-1)))
    missed_electrodes=find(brokenLeads);
    index_after_del(missed_electrodes)=[];
    [int, keepindex, repindex] = mesh_laplacian_interp(lap, index_after_del);
    data_aft_del(missed_electrodes,:)=[];
    SigInt=int*data_aft_del;
    Sig=SigInt;
    save(strcat('Run0053-ts_',num2str(aa),'_',num2str(ii),'_2'),'Sig')
end
end
