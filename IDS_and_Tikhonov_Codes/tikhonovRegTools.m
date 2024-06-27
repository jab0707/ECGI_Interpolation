function [X1,X2,lambda,singLambda] = tikhonovRegTools(Y, A, frRange)

% USAGE: [X1,X2,lambda,singLambda] = tikhonovRegTools(Y, A, frRange)
%
% Runs Tikhonov regularization for variable lambda at all time instants OR
% a single (mean or median) lambda value at all times.
% 
% INPUTS:
%
%   Y:   BSP data matrix (nLeadsTorso x nFrames)
%   A:   The forward matrix (nLeadsTorso x nLeadsHeart)
%   frRange = [init fin] is the time range for Y for which we take the median of lambda
%             (Optional)
%
% OUTPUTS:
%
%   X1: Inverse solution where a different lambda (determined by the
%   L-curve) is used at each time
%   X2: Inverse solution with a median (single) lambda value
%   lambda: The lambda values at each time
%   singLambda: median value of all lambdas within the frame range
%
% Uses Per C. Hansen's Regularization toolbox (v4)
%
% Yesim Serinagaoglu (yserin@metu.edu.tr)
% 25/06/2024

nFrames = size(Y,2);

if nargin < 3
    frRange = [1 nFrames];
end

lambda = zeros(nFrames, 1);
[X1,X2] = deal(zeros(size(A,2), nFrames));

[U,s,V] = csvd(A);

for fr = 1:nFrames
    [lambda(fr),~,~,~] = l_curve(U, s, Y(:,fr));
    if lambda(fr)>2
        if fr == 1
            lambda(fr) = 0.05;
        else
            lambda(fr)=lambda(fr-1);
        end
    end
    if lambda(fr)<0.000005
        if fr == 1
            lambda(fr) = 0.05;
        else
            lambda(fr)=lambda(fr-1);
        end
    end
    X1(:,fr) = tikhonov(U,s,V,Y(:,fr),lambda(fr));
end

singLambda = median(lambda(frRange(1):frRange(2)));
% singLambda = mean(lambda(frRange(1):frRange(2)));

for fr = 1:nFrames
    X2(:,fr) = tikhonov(U,s,V,Y(:,fr),singLambda);
end
