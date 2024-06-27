function interp_Signal = krigingInterpolation(W,knownLocations,unknownLocations,Signals)
%krigingInterpolation.m Given kriging wieghts, interpolate signals from
%electrodes to nodes


%Solving for weights - matV*W=matVS
mesh_size = size(unknownLocations,1);
for i =1:size(unknownLocations,1) 
    nn=0;
    for j = 1:size(knownLocations,1);
        nn=nn+ W(i,j)*Signals(j,:);       
    end
    interp_Signal(i,:)=nn;
end

return

