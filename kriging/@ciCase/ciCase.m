classdef ciCase
%CICASE Class to keep info on each ci case
    
%% Properties
    properties
        triTorso; % surface mesh of the electrodes on the body surface (e.g. from the vest)
        triHeart; % surface mesh of the heart (from segmentation)
        triPeri; % Specific to Mark's modeling data: exterior, fake pericardium
        normHeart; % unit vectors normal to the heart surface
        normTorso; % unit vectors normal to the torso surface
        neighsHeart; % neighbours of each vertex of the mesh
        neighsTorso; % neighbours of each vertex of the mesh
        patientData = []; % Data like, Name, link to datafiles, etc
        % Specific to Mark's modeling data, contains vm, ve, AT on triHeart, ve on
        % triTorso, ve on triPeri
        modelData = struct('vmHeart',[],'veHeart',[],'veTorso',[],'vePeri',[],'atHeart',[]); 
        landmarks = []; % Optional supplemental meshes (aorta, other arteries...)
        channelIndex; % Index of the electrodes removed (or remaining ?)
        ivHeart; % First 2 vectors of an orthonormal basis at each vertex
                 % (with the normal)
        jvHeart; % First 2 vectors of an orthonormal basis at each vertex
                 % (with the normal)
        lapHeart; % Surface Laplacian
        info; % field from Mark's model, that contains info on the experiment
    end
    properties (Dependent)
        unit;
    end
    properties (Access = private)
        punit;
    end
    methods
        %% Constructor function
        function obj = ciCase(varargin)
            function path = pathGet()
                [file, folder] = uigetfile('*.mat', 'Please select ciCase file');
                path = fullfile(folder, file);
            end
            % Parse the arguments given to ciCase constructor
            opts = parseArgsToStruct(varargin, {'path', 'heart', 'torso', 'unit','landmarks','patientData','channelIndex','torsonorm'}, [], ...
                {[], [], [], 0.01,[],[],1:252,[]});
            % Try to load the path
            if isempty(opts.heart) || isempty(opts.torso)
                % Open a gui to browse files if no file is given
                if isempty(opts.path), opts.path = pathGet(); end
                % Load the mat file
                tmp = load(opts.path, '-mat');
                % Test if the loaded structure contains .sim and .info -->
                % Mark's modeling data, otherwise --> CI data
                try
                    % Use input geometries
                    obj.triHeart = TriRep(opts.heart.Triangulation, opts.heart.X);
                    obj.triTorso = TriRep(opts.torso.Triangulation, opts.torso.X);
                     obj.channelIndex = opts.channelIndex;
                    try obj.landmarks = convertStructToTriRep(opts.landmarks); catch, end
                    if isempty(opts.patientData), obj.patientData = opts.patientData;end
                catch
                    error('Invalid input geometries');
                end
            end

            % Compute neighbors and normal vectors
            obj.neighsHeart = neighborsJD(obj.triHeart);
            obj.neighsTorso = neighborsJD(obj.triTorso);
            [obj.normHeart, obj.ivHeart, obj.jvHeart] = normalVectorsJD(obj.triHeart, [], obj.neighsHeart);
            
            if isempty(opts.torsonorm) 
                obj.normTorso = normalVectorsJD(obj.triTorso);
            else
                obj.normTorso = opts.torsonorm;
            end
            
         end

    end
end

