function obj = ciTransfer(ciCase, channelMask)
    
    assert(nargin > 0 && isa(ciCase, 'ciCase'), 'Invalid input ciCase');
    if nargin < 2 || isempty(channelMask), channelMask = true(size(ciCase.triTorso.X, 1), 1); end

    assert(any(strcmpi(formulation, obj.vformulation)), 'Invalid input formulation');
    channelMask = logical(channelMask(:));
    assert(length(channelMask) == size(ciCase.triTorso.X, 1), 'Invalid channel mask size');

    obj.pciCase = ciCase;
    obj.pchannelMask = channelMask;

    % Compute the transfer matrix
    fprintf('>> Computing transfer matrix using %s formulation...\n', formulation);
    [obj.transfer, obj.virtual, obj.virtualG, obj.virtual3DG, obj.virtualL, obj.virtual3DL] = obj.compute;
    [obj.matU,obj.vectS,obj.matV] = csvd(obj.transfer);
end