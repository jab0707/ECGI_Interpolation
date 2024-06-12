function [Beat] = BeatSegmetation(OriginalSignal,BeatsSegments)
Beat=cell(1,length(BeatsSegments));
for i=1:length(BeatsSegments)
Beat{1,i}=OriginalSignal(:,BeatsSegments{1,i}(1):BeatsSegments{1,i}(2));
end
end

