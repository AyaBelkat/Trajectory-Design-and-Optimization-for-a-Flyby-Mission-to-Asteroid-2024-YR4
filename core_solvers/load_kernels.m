function load_kernels()
    % load_kernels: Loads SPICE ephemeris and data for the 2024 YR4 mission.
    % This utility ensures all required kernels are loaded from the /data folder
    % using relative path resolution via the 'which' command.

    % 1. Clear existing kernels to prevent data contamination
    cspice_kclear;

    % 2. Define the list of kernels and meta-kernels to load
    % (Based on your verified data list)
    kernelsToLoad = { ...
        'naif0012.tls.pc', ...            % Leap seconds
        'pck0011.tpc', ...                % Planetary constants
        '54509621_2024yr4_data.bsp', ...  % Custom Asteroid SPK
        'mission_meta.tm', ...            % Primary Meta-Kernel
        'mission_equa_data.tm'            % Secondary Data Map
    };

    fprintf('Loading SPICE Kernels...\n');

    % 3. Loop through and furnish each file
    for i = 1:length(kernelsToLoad)
        currentFile = kernelsToLoad{i};
        fullPath = which(currentFile); % Resolves the path inside /data

        if ~isempty(fullPath)
            cspice_furnsh(fullPath);
            fprintf('  [OK]  %s\n', currentFile);
        else
            % If it's the large ephemeris file, provide a specific warning
            if strcmp(currentFile, 'de430.bsp')
                warning('Large ephemeris file (de430.bsp) missing from /data folder.');
            else
                warning('Kernel not found: %s. Check /data folder.', currentFile);
            end
        end
    end
    
    % Note: Large files like de430.bsp should be loaded via mission_meta.tm 
    % or explicitly added here if the user has downloaded them.
    fprintf('SPICE Environment Ready.\n');
end
