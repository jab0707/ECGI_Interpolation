function [ dSi, dSj, iVect, jVect ] = dVdX( obj, sig )
%dVdX approximates spatial derivative of the potentials along the geometry
% Input: 
%       obj ciCase
%       sig : input signal
%
% Output:
%       dSi = gradient along direction i
%       dSj = gradient along direction j
%       iVect = unit vectors i coordinates
%       jVect = unit vectors j coordinates
%
% JD 12/11/14, modified 23/06/15

%% Input management
    [nPts, nSamp] = size(sig);
    assert(nPts == size(obj.triHeart.X, 1), 'Invalid signal dimensions')
    
%% Calculate unit vectors from first neighbor and normal vectors
    iVect = nan(nPts, 3);
    jVect = nan(nPts, 3);
    goodPts = ~cellfun(@isempty, obj.neighsHeart)';
    for i = find(goodPts)
        vect1 = obj.triHeart.X(obj.neighsHeart{i}(1),:) - obj.triHeart.X(i,1);
        jVect(i,:) = cross(vect1, obj.normHeart(i,:));
        jVect(i,:) = jVect(i,:) ./ norm(jVect(i,:));
        iVect(i,:) = cross(obj.normHeart(i,:), jVect(i,:));
    end
    
%% Calculate gradients for each point
    dSi = zeros(nPts, nSamp);
    dSj = zeros(nPts, nSamp);
    for i = find(goodPts)
        nNeighs = length(obj.neighsHeart{i});
        % Calculate signal difference with neighbors
        dSig = sig(obj.neighsHeart{i},:) - repmat(sig(i,:), nNeighs, 1);
        
        % Calculate vectors from point to neighbors
        dVect = obj.triHeart.X(obj.neighsHeart{i},:) - repmat(obj.triHeart.X(i,:), nNeighs, 1);
        
        % Calculate the projection of this difference on the unit vectors
        iComp = dot(dVect, repmat(iVect(i,:), nNeighs, 1), 2); % i components
        jComp = dot(dVect, repmat(iVect(i,:), nNeighs, 1), 2); % j components
        
        diSig = dSig .* repmat(iComp, 1, nSamp);
        diSig = diSig ./ repmat(sqrt(sum(dVect.^2, 2)), 1, nSamp);
        djSig = dSig .* repmat(jComp, 1, nSamp);
        djSig = djSig ./ repmat(sqrt(sum(dVect.^2, 2)), 1, nSamp);
        
        % Solve in least squares sense the components
        dSi(i,:) = ones(nNeighs,1)\diSig;
        dSj(i,:) = ones(nNeighs,1)\djSig;
    end
end


