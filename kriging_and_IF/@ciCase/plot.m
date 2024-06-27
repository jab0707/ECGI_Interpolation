function hax = plot(obj, data, varargin)
% Modified by JD 06/01/15

    if nargin < 2, data = []; end
    if ischar(data), 
        varargin = cat(2, {data}, varargin);
        data = [];
    end
    opts = parseArgsToStruct(varargin, {'sat', 'cmap', 'fitmap', 'dVdT', 'mapItv'}, [], {[], 'jet', 1, -1, [1 size(data, 2)]});
    
    %% INPUT MANAGEMENT
    if numel(data) ~= length(data)
        assert(size(data, 1) == size(obj.triHeart.X, 1), 'Data & mesh mismatch');
        switch(opts.dVdT)
            case 0
                error('Invalid dVdT argument');
            case 1 % Use max dVdT
                data = dVdT(data, 0, [], opts.mapItv);
            case -1 % Use min dVdT
                data = dVdT(data, 1, [], opts.mapItv);
        end
    end
    
    data = double(data(:));
    s = opts.sat(:);
    assert(isempty(s) || all(size(s) == size(data)), 'Color & sat dimensions mismatch');
    assert(isempty(data) || size(obj.triHeart.X, 1) == length(data), 'Data & mesh dimensions mismatch');
    
    %% Calculate color
    if isempty(data)
        c = obj.triHeart.X(:, 3);
        cmfun = str2func(opts.cmap);
    elseif isempty(s)
        c = data;
        cmfun = str2func(opts.cmap);
        cmfun = @(x) flipud(cmfun(x));
    else
        h = - double(data);
        h = h - min(h);
        h = 127 * h / max(h);
        h = round(h+1);
        
        cmfun = str2func(opts.cmap);
        color = cmfun(128);
        cmfun = @(x) flipud(cmfun(x));
        
        h = color(h, :);
        
        if min(s) < 0, s = s - min(s); end
        if max(s) > 0, s = s / max(s); end
        
        % Fade to gray in case of low input saturation
        c = h .* repmat(s, 1, 3) + repmat(0.5-s/2, 1, 3);
    end
    
    hax = trisurf(obj.triHeart, 'EdgeAlpha',0,'SpecularColorReflectance',0,'AmbientStrength',0.4,'FaceColor','interp', 'FaceVertexCData', c);
    colormap(cmfun(128));
    if opts.fitmap && ~isempty(data), caxis([min(data), max(data)]); end
    if ~isempty(data), colorbar(); end

    % Plot the landmarks
    if ~isempty(obj.landmarks)
        lnames = fieldnames(obj.landmarks);
        hold on
        for i=1:length(lnames)
            trisurf(obj.landmarks.(lnames{i}), 'EdgeAlpha',0, 'FaceColor', [0 0 0], 'FaceAlpha', 0.2);
        end
        hold off
    end
    
    axis off
    axis equal
    axis tight
end

