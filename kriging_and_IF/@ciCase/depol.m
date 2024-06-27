function [ map ] = depol( obj, sig, varargin )
%DEPOL Compute depolarization map using different techniques

    %% Parse input
    % Varargin
    opts = parseArgsToStruct(varargin, ...
        {'method', 'thresh1', 'thresh2', 'remap', 'delayfun', 'mapItv'}, [], ...
        {'weighted', [], [], 0, 'hcorrd', [1 size(sig, 2)]});
    
    % Check input signal size
    nPts = size(sig, 1);
    assert(nPts == size(obj.triHeart.X, 1), 'Invalid input signal dimensions');
    sig = sig(:, opts.mapItv(1):opts.mapItv(2));
    
    %% Compute dVdT and confidence values
    [gdVdT, singleConf] = arrayfun(@(x) uniLAT2(sig(x,:)), (1:nPts)');
    
    %% Compute delay matrix and confidence values
    [lagMat,delayConf] = delayMat(obj, sig, 'delayfun', opts.delayfun);

    %% Build and solve equation system
    switch lower(opts.method)
        case 'thresh'
            map = thresholdEquations(gdVdT, singleConf, lagMat, delayConf, obj.triHeart, obj.neighsHeart, opts);
        case 'weighted'
            map = weightedEquations(gdVdT, singleConf, lagMat, delayConf, obj.triHeart, obj.neighsHeart, opts);
        case 'equal'
            map = equalEquations(gdVdT, singleConf, lagMat, delayConf, obj.triHeart, obj.neighsHeart, opts);
        case 'curve'
            map = curveEquations(gdVdT, singleConf, lagMat, delayConf, obj.triHeart, obj.neighsHeart, opts);
        otherwise
            error('Invalid input method "%s"', opts.method);
    end
    
    %% Remap if necessary
    if opts.remap
        map = remap(map, sig);
    end
    
    %% Correction for mapping interval
    map = map + opts.mapItv(1) - 1;
end

%% EQUAL WEIGHTS FUNCTION
function map = equalEquations(lat, ~, delay, delayConf, tri, ~, ~)
    nPts = size(tri.X, 1);
    [r, c] = find(delayConf);
    delayList = [r c];
    delayList = sort(delayList, 2);
    delayList = unique(delayList, 'rows');
    nDelay = size(delayList, 1);

    % Write sparse equation system
    eqSys = sparse([1:nDelay, 1:nDelay, (nDelay+1):(nDelay+nPts)], ...
        [delayList(:,1)', delayList(:,2)' 1:nPts], ...
        [ones(1, nDelay), -ones(1, nDelay), ones(1, nPts)], ...
        nDelay + nPts, nPts);

    % Write constants
    eqConst = ...
        [delay(sub2ind(size(delay), delayList(:, 1), delayList(:, 2))) ;...
        lat];

    % Solve using least-squares
    map = full(eqSys \ eqConst);
end

%% CURVE FUNCTION
function map = curveEquations(lat, latConf, delay, delayConf, tri, ~, ~)
    nPts = size(tri.X, 1);
    [r, c] = find(delayConf);
    delayList = [r c];
    nDelay = size(delayList, 1);

    % Write sparse equation system
    eqSys = sparse([1:nDelay, 1:nDelay, (nDelay+1):(nDelay+nPts)], ...
        [delayList(:,1)', delayList(:,2)' 1:nPts], ...
        [ones(1, nDelay), -ones(1, nDelay), ones(1, nPts)], ...
        nDelay + nPts, nPts);

    % Write constants
    eqConst = ...
        [delay(sub2ind(size(delay), delayList(:, 1), delayList(:, 2))) ;...
        lat];    
    
    latWFun = @(c) 1./exp(-98.3 .* c + 7.33); % Model adjusted parameters
    delayWFun = @(c) 1./exp(-6.05 .* c + 8.40);
    
    latW = latWFun(latConf);
    delayW = full(delayConf(sub2ind(size(delayConf), delayList(:, 1), delayList(:, 2)))); 
    delayW = delayWFun(delayW);
    
    % Build curve
    nLambda = 100;
    %lspace = logspace(-2, 2, nLambda);
    lspace = linspace(0.01, 0.99, nLambda);
    
    latResiduals = zeros(nLambda, 1);
    delayResiduals = zeros(nLambda, 1);
    
    delayOnlyMat = sparse([1:nDelay, 1:nDelay], ...
        [delayList(:,1)', delayList(:,2)'], ...
        [ones(1, nDelay), -ones(1, nDelay)], ...
        nDelay, nPts);
    delayOnly = delay(sub2ind(size(delay), delayList(:, 1), delayList(:, 2)));
    map = zeros(nPts, nLambda);
    for i = 1:nLambda
        %eqWeight = [lspace(i)*ones(nDelay, 1); ones(nPts, 1)];
        %eqWeight = [lspace(i)*ones(nDelay, 1)./nDelay; (1-lspace(i))*latW./sum(latW)];
        
        eqWeight = [lspace(i)*ones(nDelay, 1); (1-lspace(i))*ones(nPts, 1)];
        %eqWeight = [lspace(i)*delayW; (1-lspace(i))*latW];
        
        map(:,i) = lscov(eqSys, eqConst, eqWeight);
        latResiduals(i) = sqrt(mean((lat - map(:,i)).^2));
        delayResiduals(i) = sqrt(mean((delayOnlyMat*map(:,i) - delayOnly).^2));
    end
    
%     figure('Name', 'Residuals', 'Units', 'Normalized', 'Position', [0.2 0.2 0.6 0.6 ])
%     subplot(1,2,1);
%     %semilogx(lspace, [latResiduals delayResiduals latResiduals+delayResiduals ]); 
%     plot(lspace, [latResiduals delayResiduals latResiduals+delayResiduals]); 
%     legend({'RMS LAT residual', 'RMS delay residual', 'Sum of residuals'});
%     title('Residuals');
%     xlabel('\lambda');
%     
%     subplot(1,2,2);
%     %semilogx(lspace, [curvature( latResiduals, lspace) curvature(delayResiduals, lspace)]);
%     plot(lspace, [curvature( latResiduals, lspace) curvature(delayResiduals, lspace)]);
%     legend({'RMS LAT residual', 'RMS delay residual'});
%     title('Curvature of residuals');
%     xlabel('\lambda');
%     
%     figure()
%     semilogx(delayResiduals, curvature(delayResiduals, latResiduals));
    figure()
    plot(lspace, [latResiduals delayResiduals latResiduals+delayResiduals]); 
end

%% THRESHOLD FUNCTION
function map = thresholdEquations(lat, latConf, delay, delayConf, tri, neighbors, opts)
    nPts = size(tri.X, 1);
    if isempty(opts.thresh1), opts.thresh1 = 10; end
    if isempty(opts.thresh2), opts.thresh2 = 0.7; end

    % Find unipolars above threshold
    latList = latConf < opts.thresh1;
    latList = filterMesh(tri, latList, 'erode', neighbors);
    latList = find(latList);
    nLat = size(latList, 1);

    % Find bipolars above threshold
    [r, c] = find(delayConf > opts.thresh2);
    delayList = [r c];
    delayList = sort(delayList, 2);
    delayList = unique(delayList, 'rows');
    nDelay = size(delayList, 1);

    % Write sparse equation system
    eqSys = sparse([1:nDelay, 1:nDelay, (nDelay+1):(nDelay+nLat)], ...
        [delayList(:,1)', delayList(:,2)' latList'], ...
        [ones(1, nDelay), -ones(1, nDelay), ones(1, nLat)], ...
        nDelay + nLat, nPts);

    % Write constants
    eqConst = ...
        [delay(sub2ind(size(delay), delayList(:, 1), delayList(:, 2))) ;...
        lat(latList)];

    % Check for undetermined points
    badPts = checkIsolation(eqSys, latList, delayList);

    % Solve system using least-squares
    map = NaN(nPts,1);
    map(~badPts) = eqSys(:, ~badPts) \ full(eqConst);

end

%% WEIGHTED LEAST SQUARES FUNCTION
function map = weightedEquations(lat, latConf, delay, delayConf, tri, ~, opts)
    nPts = size(tri.X, 1);
    % Default (model adjusted) weighing functions
    if isempty(opts.thresh1), opts.thresh1 = @(c) 1./exp(-98.3 .* c + 7.33); end
    if isempty(opts.thresh2), opts.thresh2 = @(c) 1./exp(-6.05 .* c + 8.40); end
    
    assert(isa(opts.thresh1, 'function_handle'), 'Invalid input LAT weighing function');
    assert(isa(opts.thresh2, 'function_handle'), 'Invalid input delay weighing function');
    
    latWFun = opts.thresh1;
    delayWFun = opts.thresh2;

    % Find delays 
    [r, c] = find(delayConf);
    delayList = [r c];
    nDelay = size(delayList, 1);

    % Write sparse equation system
    eqSys = sparse([1:nDelay, 1:nDelay, (nDelay+1):(nDelay+nPts)], ...
        [delayList(:,1)', delayList(:,2)' 1:nPts], ...
        [ones(1, nDelay), -ones(1, nDelay), ones(1, nPts)], ...
        nDelay + nPts, nPts);

    % Write constants
    eqConst = ...
        [delay(sub2ind(size(delay), delayList(:, 1), delayList(:, 2))) ;...
        lat];

    % Write equation weights
    latW = latWFun(latConf);
    delayW = full(delayConf(sub2ind(size(delayConf), delayList(:, 1), delayList(:, 2))));
    delayW = delayWFun(delayW);
    eqWeight = [delayW ; latW];

    % Solve system using weighted least squares
    map = lscov(eqSys, eqConst, eqWeight);
end

%% REMAPPING FUNCTION
function map = remap(map, sig)
    % Find candidates for remaping
    cands = diff(sig, [], 2) < 0.3 * repmat(min(diff(sig, [], 2), [], 2), 1, size(sig, 2) - 1);

    % Look for candidate closest to eqSol
    for i = 1:size(sig, 1)
        try
            icands = find(cands(i,:));
            [~, ibest] = min(abs(icands - map(i)));
            map(i) = icands(ibest);
        catch
            fprintf('Could not remap pt n°%d\n', i);
        end
    end
end

%% CHECK ISOLATION FUNCTION
function badPts = checkIsolation(inSys, singleList, delayList)
    % First find points completely isolated
    badPts = ~sum(abs(inSys));
%     if any(badPts)
%         fprintf('Found %d isolated (undetermined) point(s)\n', sum(badPts));
%     end
    
    % Next go through data and region grow connected points (by delay)
    nRegions = 0;
    regions = {};
    while 1
        if isempty(delayList), break; end
        nRegions = nRegions+1;
        [regions{nRegions}, delayList] = growDelays(delayList(1,:)', delayList(2:end, :));
    end
    
%     if length(regions) > 1
%         fprintf('%d seperate regions\n', length(regions));
%     end
    
    % Verify that all regions have activation times
    for i=1:length(regions)
        if isempty(intersect(regions{i}, singleList)) 
%             fprintf('Found a disconnected region of %d points with no ref. AT\n', ...
%                 length(regions{i}));
            badPts(regions{i}) = true;
        end
    end
end

%% GROWDELAYS FUNCTION
function [okPts, delayList] = growDelays(okPts, delayList)
    % Find connected points
    [~, ia, ~] = intersect(delayList(:, 1), okPts);
    [~, ib, ~] = intersect(delayList(:, 2), okPts);
    ia = unique([ia;ib]);
    if ~isempty(ia)
        cand = delayList(ia, :);
        delayList(ia,:) = [];
        okPts = unique([okPts;cand(:)]);
        
        % Recursive call to grow function
        [okPts, delayList] = growDelays(okPts, delayList);
    end
end