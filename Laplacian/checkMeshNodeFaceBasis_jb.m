function [mesh,nodeNames,elemNames] = checkMeshNodeFaceBasis_jb(mesh,~,flag,nodeName,elemName)
    %% HELP:
    % Checks the input struct mesh to ensure that the node and element fields
    % are in either row or column major orientation based on the input flag.
    % nodes can be in fields named   : 'node','pts','pt','vertices','nodes','points'
    % elements can be in fields named: 'elem','element','cell','tri','face','faces','cells','fac'
    % or the node name and elem name can be provided.
    % author: Jake Bergquist <jake.a.bergquist@gmail.com>
    meshFields = fieldnames(mesh);
    if ~exist('flag','var');flag = 'row';end

    if ~exist('nodeName','var')
        nodeNames = {};
        for nodeTest = {'node','pts','pt','vertices','nodes','points'}
            if any(strcmp(meshFields,nodeTest{1}))
                nodeNames = [nodeNames,nodeTest];
            end
        end
    else
        nodeNames = {nodeName};
    end
    if ~exist('elemName','var')
        elemNames = {};
        for elemTest = {'elem','element','cell','tri','face','faces','cells','fac'}
            if any(strcmp(meshFields,elemTest{1}))
                elemNames = [elemNames,elemTest];
            end
        end
    else
        elemNames = {elemName};
    end
    

    switch flag
        case 'row'
            dimCheck = 1;
            otherDim = 2;
        case 'col'
            dimCheck = 2;
            otherDim = 1;
    end
    
    for elemCheck = elemNames
        elemCheck = elemCheck{1};
        mesh = checkElems(mesh,elemCheck,dimCheck,otherDim);
    end
    for nodeCheck = nodeNames
        nodeCheck = nodeCheck{1};
        mesh = checkNodes(mesh,nodeCheck,dimCheck,otherDim);
    end

end

function mesh = checkElems(mesh,elemName,dimCheck,otherDim)
    
    elemSize = size(mesh.(elemName));
    if max(elemSize) <= 4
        warning('Meshes with 4 or fewer elements may not be transposed correctly\n')
    end
    elemNumber = min(elemSize);
    if ~strcmp(elemName,'')
        if elemSize(dimCheck) == elemNumber
            %make sure it is an int64
            mesh.(elemName) = int64(mesh.(elemName));
        elseif elemSize(otherDim) == elemNumber
            mesh.(elemName) = int64(mesh.(elemName))';
        else
            error('Dim check failed, somethign wrong with the mesh nodes. Elems not in expected format. Expected %d x n, got %d x %d',elemNumber,elemSize(1),elemSize(2))
        end
    end

end


function mesh = checkNodes(mesh,nodeName,dimCheck,otherDim)
    nodeNumber = 3;
    nodeSize = size(mesh.(nodeName));
    if ~strcmp(nodeName,'')
        if nodeSize(dimCheck) == nodeNumber
            %do nothing
        elseif nodeSize(otherDim) == nodeNumber
            mesh.(nodeName) = mesh.(nodeName)';
        else
            error('Dim check failed, somethign wrong with the mesh nodes. Nodes not in expected format. Expected %d x n, got %d x %d',nodeNumber,nodeSize(1),nodeSize(2))
        end
    end

end