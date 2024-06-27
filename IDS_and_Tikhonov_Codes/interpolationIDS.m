function [complete_data] = interpolationIDS(reduced_data, geompts, brokenLeads, R)

% USAGE: [epidata] = interp_epi(sockdata,sockgeom,epigeom);
%
% INPUTS:
%
%   reduced_data:   Data matrix including only the good leads (nLeadsReduced x nFrames)
%   geompts:        Sock node points (nLeads x 3)
%   brokenLeads:    A vector of size nLeads x 1. Equals 1 for broken leads, 0 for the good leads.  
%   R:              The radius around the bad lead. If provided, use the
%                   good leads inside this radius for interpolation. If not given, all
%                   good leads are used.
%
% OUTPUT:
%
%   complete_data:        A matrix of interpolated data (nLeads x nFrames)
%
% Uses inverse-distance-squared interpolation:
%
% D. Shepard, ‘A two-dimensional interpolation function for
% irregularly-spaced data,’ Proceedings of ACM National Conference, 1968.
% 
% Yesim Serinagaoglu (yserin@metu.edu.tr)
% 25/06/2024

badLeadIndices = find(brokenLeads==1);
goodLeadIndices = find(brokenLeads==0);

ptsGoodLeads = geompts(goodLeadIndices,:);
ptsBadLeads = geompts(badLeadIndices,:);

nGoodLeads = length(goodLeadIndices);
nBadLeads = length(badLeadIndices);

% Get sizes of data
% nLeads = size(geompts,1);
nFrames = size(reduced_data,2);

% Initializations
distance = zeros(nBadLeads,nGoodLeads);

% Apply interpolation to obtain EpiData
interpolated = zeros(nBadLeads,nFrames);
for i = 1:nBadLeads
    for j = 1:nGoodLeads
        distance(i,j) = sqrt((ptsBadLeads(i,1)-ptsGoodLeads(j,1))^2 + (ptsBadLeads(i,2)-ptsGoodLeads(j,2))^2 + (ptsBadLeads(i,3)-ptsGoodLeads(j,3))^2);
%         inv_distance2(i,j) = 1/(distance(i,j)^2);
    end

    if nargin == 4 % R is provided
        % If you want to include only those leads within a radius
        includeLeads = (find(distance(i,:) <= R));
        while length(includeLeads)<=6 % make sure there are at least 6 good leads for interpolation
            R=R+5;
            includeLeads = (find(distance(i,:) <= R));
        end
    
        newdistance=distance(i,includeLeads);
        inv_distance2 = 1./(newdistance.^2);
        num_sum_vec = inv_distance2*reduced_data(includeLeads,:);
        denom_sum = sum(inv_distance2); % row-wise summation
    else
        % If you want to include all leads
        inv_distance2 = 1./(distance(i,:).^2);
        num_sum_vec = inv_distance2*reduced_data;
        denom_sum = sum(inv_distance2); % row-wise summation
    end

    interpolated(i,:) = num_sum_vec/denom_sum;
end

complete_data(badLeadIndices,:)=interpolated;
complete_data(goodLeadIndices,:)=reduced_data;
