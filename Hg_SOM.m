% SOM Tool Box original files at http://www.cis.hut.fi/somtoolbox/download/, click on [SOM Toolbox 2.0 (Mar 17 2005)]
% The "somtoolbox" folder is identical to the original except for the file "som_make" which now allows algorithm, initalg, mapsize and training to be selected from here
% The experiments were performed using "Octave-10.3.0" during the October of 2025.
% List of the files:
%   Cossarini_Model_Output/
%             20200101_y-CMCC--PSAL-MFSe3r1-MED-b20220901_re-sv01.00.nc
%             20200101_y-CMCC--TEMP-MFSe3r1-MED-b20220901_re-sv01.00.nc
%             20200101_y-OGS--BIOL-MedBFM3-MED-b20221120_re-sv05.00.nc
%   Rosati_Data_Model_Output/
%             ave.20140616-00_00_00.DMHg.nc
%             ave.20140616-00_00_00.MMHg.nc
%             data_ancillary_DGM.csv
%             data_ancillary_HgT.csv
%             data_ancillary_MeHg.csv
%   somtoolbox/
%             SOM original scripts
%   Puglia_Hg_SOM.m (THIS FILE)

%% ============================================================
%% ========================= START ============================
%% ============================================================
clc                                                                                                  % Clean command window
clear                                                                                                % Clean workspace
DIR = "C:/Users/mpuglia/OneDrive - OGS/Desktop/Puglia_Bandelj_Cossarini_Rosati_Experiment_10-2025";  % Main Path
cd(DIR)                                                                                              % Set directory
addpath('somtoolbox');                                                                               % Path to the original functions




%% ============================================================
%% ================= OCTAVE Needed Packages ===================
%% ============================================================
%pkg install -forge netcdf         % Install the NetCDF package (uncomment if not installed)
pkg load netcdf                    % Load the NetCDF package in Octave
pkg load io                        % Only needed if the io package is not already loaded




%% ============================================================
%% ================ Rosati's Hg Observations ==================
%% ============================================================

%% ===================== 0. Load data =====================
dataset = csvread('Rosati_Data_Model_Output/data_ancillary_MeHg.csv', 1, 0);     % Load dataset
datasetval = dataset(:, [9 10 11 12]);                                           % Select MeHg/HgT, Temp, Sal, O2 observations
%datasetval = dataset(:, [11 16 17 18]);                                         % Select DGM, Temp, Sal, O2 observations

%% ===================== 1. Data into structure =====================
sData = som_data_struct(datasetval, 'name', 'MeHg Dataset', 'comp_names', {'MeHg','Temp','Sal','O2_uM'});  % Wrap in SOM Toolbox data structure


%% ===================== 2. Normalize data =====================

sData = som_normalize(sData, 'var');                    % Normalize data


%% ===================== 3. Create, initialize and train a SOM. =====================
% 'algorithm'  *(string) training: 'seq' or 'batch' (default) or 'sompak'
% 'init'       *(string) initialization: 'randinit' or 'lininit' (default)
% 'mapsize'    *(string) do you want a 'small', 'normal' or 'big' map
% 'training'    (string) 'short', 'default' or 'long'
algorithm = 'batch';                                              % Selected: 'seq'/'batch'
initalg = 'lininit';                                            % Selected: 'lininit'
mapsize = 'big';                                                % Selected: 'big'
training = 'long';                                              % Selected: 'long'
sMap = som_make(sData, algorithm, initalg, mapsize, training);  % Train the SOM with sData and specified options

clear algorithm initalg mapsize training                        % Clear temporary variables from workspace

%% ===================== 4. Find Final BEST-MATCHING UNITS (BMUs) =====================
bmus = som_bmus(sMap, sData);                           % Find the Best-Matching Unit (BMU) for each data sample in sData using the trained SOM


%% ===================== 5. SOM quality metrics =====================
[q, t] = som_quality(sMap, sData);                      % Compute SOM quality metrics: q = quantization error, t = topographic error


%% ===================== 6. Save results =====================
save('C:\Users\mpuglia\MeHg_Rosati_SOM_Observations.mat')


%% =====================================================
%% ===== PLOTS ===== PLOTS ===== PLOTS ===== PLOTS =====
%% =====================================================

%% ===================== 7. Visualize component planes =====================
figure(1);                                       % Open figure

colormap(1-gray)                                 % Use an inverted grayscale colormap for better contrast
som_show(sMap, 'norm', 'd');                     % Display the SOM component planes normalized by distance
text(1, -7, 'Component planes for MeHg', ...     % Add a title text to the figure at position (1, -7)
     'HorizontalAlignment', 'center', ...        % Center the text horizontally
     'FontSize', 14, 'FontWeight', 'bold');      % Use larger, bold font for emphasis
saveas(gcf, 'MeHg_SOM_component_planes.png');    % Save as PNG


%% ===================== 8. Hit Plots =====================
hit = som_hits(sMap, sData);                     % Compute the number of data samples (hits) mapped to each SOM node

% ====== Plot marks based on size ======
figure(2);                                             % Open figure

colormap(1-gray)                                       % Set colormap to inverted grayscale for better contrast
som_show(sMap, 'norm', 'd', 'umat', 'all');            % Display only the U-Matrix (distance map) of the SOM
som_show_add('hit', hit, ...                           % Overlay hit markers on the SOM map
             'MarkerColor', 'r', 'MarkerSize', 0.5);   % Set marker color to red and marker size to 0.5
title('SOM Hits Sized Marks');                         % Add a title to the plot
saveas(gcf, 'MeHg_SOM_hit_mark.png');       % Save as PNG

% ====== Plot numeric hit counts ======
figure(3);                                                                   % Open figure

colormap(1-gray)                                                             % Set colormap to inverted grayscale for better contrast
som_show(sMap, 'norm', 'd', 'umat', 'all');                                  % Display only the U-Matrix (distance map) of the SOM
pos = som_unit_coords(sMap);                                                 % Get the (x, y) coordinates of all SOM units (neurons)

x_offset = 0.9;                                                              % Horizontal offset
y_offset = 0.9;                                                              % Vertical offset
pos(:,1) = pos(:,1) * sqrt(3)/1.73;                                          % Horizontal scale
pos(:,2) = pos(:,2) * 1.155;                                                 % Vertical scale

for i = 1:length(hit)                                                        % Loop through each SOM unit
    if hit(i) > 0                                                            % Only label units that have one or more hits
        text(pos(i,1) + x_offset, pos(i,2) + y_offset, num2str(hit(i)), ...  % Place text near the unit showing the hit count
              'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');           % Set text color, size, and bold font for visibility
    end
end
title('SOM Hit Numbers');                                                    % Add a title to the plot
saveas(gcf, 'MeHg_SOM_hit_number.png');       % Save as PNG


% ====== Plot regional labels ======
figure(4);                                                                   % Open figure

textdata = csv2cell('Rosati_Data_Model_Output/data_ancillary_MeHg.csv');     % Read CSV file into a cell array
RegionalLabels = textdata(2:end,5);                                          % Extract the regional labels from the 5th column (skip header row)
sData.labels = RegionalLabels;                                               % Initialize your SOM data structure labels with regional labels
map_keys   = {'SAD','TYR2','ION1_2','LEV','SWM','ALB','NWM','ION3','TYR1'};  % Original CSV labels
map_values = {  'D',   'T',     'I',  'L',  'S',  'A',  'N',   'I',   'T'};  % Corresponding single-letter codes

#RegionalLabels = textdata(2:end,16);                                         % Extract the seasonal labels from the 5th column (skip header row)
#sData.labels = RegionalLabels;                                               % Initialize your SOM data structure labels with seasonal labels
#map_keys   = {'autumn','spring','summer'};                                   % Original CSV labels
#map_values = {  'A',   'P',     'U'};                                        % Corresponding single-letter codes

mappedLabels = RegionalLabels;                                               % Initialize mappedLabels with the original labels

% Loop through each key in map_keys
for k = 1:length(map_keys)
    idx = strcmp(RegionalLabels, map_keys{k});                               % Find all occurrences of the current key in RegionalLabels
    mappedLabels(idx) = map_values(k);                                       % Replace matching labels with the corresponding single-letter code
end

sData.labels = mappedLabels;                                                 % Assign the mapped labels to the SOM data structure
sMap = som_autolabel(sMap, sData, 'vote');                                   % Automatically label the SOM map using the mapped labels
colormap(1-gray)                                                             % Set colormap to grayscale (1-gray inverts the grayscale)
som_show(sMap,'norm','d', 'umat', 'all');                                    % Display the SOM map with normalized distances
som_show_add('label', sMap, 'TextSize', 10, 'TextColor', 'r');               % Add labels normally
hText = findall(gca, 'Type', 'text');                                        % Get all text objects in the current axes
set(hText, 'FontWeight', 'bold');                                            % Set font weight to bold
title('SOM Regional Labels');                                                % Add a title to the plot
saveas(gcf, 'MeHg_SOM_hit_region.png');       % Save as PNG


%% ===================== 9. Map SOM units back to coordinates =====================    % NOT USEFUL OUTPUT PRODUCED
figure(5);                                                                            % Open figure

top = 100                                                                             % Select the top of the water column where to find the sample
bottom = 120                                                                          % Select the bottom of the water column where to find the sample
surface_idx = dataset(:,8) >= top & dataset(:,8) <= bottom;                           % Logical filter index for depths ># and <=#
coords_surface = dataset(surface_idx, [4 3]);                                         % Extract Longitude (column 4) and Latitude (column 3) for surface points
bmus_surface =  dataset(surface_idx,9);                                                     % BMU values corresponding to surface points

scatter(coords_surface(:,1), coords_surface(:,2), 150, bmus_surface, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);        % Scatter plot of surface points colored by BMU; 50 = marker size, 'filled' = filled circles
axis equal;                                                                           % Ensure x and y axes have the same scale
xlim([-5 30]);                                                                        % Fixed longitude range (x-axis) from -5 to 30
ylim([30 45]);                                                                        % Fixed latitude range (y-axis) from 35 to 45
colormap(viridis);
caxis([0 0.4]);
colorbar;                                                                             % Show colorbar representing BMU indices
xlabel('Longitude');                                                                  % Label x-axis
ylabel('Latitude');                                                                   % Label y-axis
title(sprintf('SOM Cluster Mediterranean Map - Depth %.1f to %.1f m', top, bottom));  % Add a title to the plot
set(gcf, 'Position', [100, 100, 800, 250]);                                           % Resize figure window for better clarity

clear map_keys map_values mappedLabels x_offset y_offset     % Clear temporary variables from workspace




%% ============================================================
%% ================== Cossarini's Model Data ==================
%% ============================================================
idx_list = [1, 5, 8, 11, 14, 16, 26]      % Loop over the slice indexes at wich you want to make the plot

%% ===================== 1. Initialise =====================
baseFolder = 'Cossarini_Model_Output';                                             % Base folder containing NetCDF files
fileNames = { ...                                                                  % List of NetCDF file names to process
    '20200101_y-CMCC--PSAL-MFSe3r1-MED-b20220901_re-sv01.00.nc', ...
    '20200101_y-CMCC--TEMP-MFSe3r1-MED-b20220901_re-sv01.00.nc', ...
    '20200101_y-OGS--BIOL-MedBFM3-MED-b20221120_re-sv05.00.nc' ...
};

files = cellfun(@(f) fullfile(baseFolder, f), fileNames, 'UniformOutput', false);  % Create full paths for all files
allVarNames = cell(1, numel(files));                                               % Initialize cell array to store variable names for each file

for f = 1:numel(files)                                                             % Loop over each NetCDF file
    ncid = netcdf.open(files{f}, 'NC_NOWRITE');                                    % Open file in read-only mode
    [~, nvars, ~, ~] = netcdf.inq(ncid);                                           % Get the number of variables in the file
    varnames = cell(1, nvars);                                                     % Initialize cell array for variable names
    for i = 0:nvars-1                                                              % Loop through all variables
        [varname, ~, ~, ~] = netcdf.inqVar(ncid, i);                               % Get the variable name
        varnames{i+1} = varname;                                                   % Store in the cell array (MATLAB/Octave indexing starts at 1)
    end

    allVarNames{f} = varnames;                                                     % Store variable names for this file
    netcdf.close(ncid);                                                            % Close the NetCDF file
end

% Display variable names for each file
for f = 1:numel(files)                                                             % Loop over each file
    fprintf('Variables in file %d:\n', f);                                         % Print file index
    disp(allVarNames{f});                                                          % Display the variable names
end

clear f i ncid nvars varname varnames                                              % Clear temporary variables from workspace


%% ===================== 2. Load .nc files =====================
ncid = netcdf.open(files{1}, 'NC_NOWRITE');      % Open first NetCDF file (PSAL) in read-only mode
varid = netcdf.inqVarID(ncid, 'so');             % Get variable ID for 'so' (salinity)
dataso = netcdf.getVar(ncid, varid);             % Read the 'so' variable data into MATLAB/Octave
netcdf.close(ncid);                              % Close the first NetCDF file

ncid = netcdf.open(files{2}, 'NC_NOWRITE');      % Open second NetCDF file (TEMP) in read-only mode
varid = netcdf.inqVarID(ncid, 'thetao');         % Get variable ID for 'thetao' (temperature)
dataT = netcdf.getVar(ncid, varid);              % Read the 'thetao' variable data
netcdf.close(ncid);                              % Close the second NetCDF file

ncid = netcdf.open(files{3}, 'NC_NOWRITE');      % Open third NetCDF file (BIOL) in read-only mode
varid = netcdf.inqVarID(ncid, 'o2');             % Get variable ID for 'o2' (oxygen)
dataO = netcdf.getVar(ncid, varid);              % Read the 'o2' variable data
netcdf.close(ncid);                              % Close the third NetCDF file

clear ncid varid                                 % Clear temporary NetCDF handles from workspace

dataT = dataT(12:end, :, 1:125);                                    % Subset temperature data: skip first 11 rows, keep all columns, limit to first 125 depth levels
dataso = dataso(12:end, :, 1:125);                                  % Subset salinity data similarly


for idx = idx_list
    %% ===================== 3. Process the matrix =====================

    slice   = [1, 5, 8, 11, 14, 16, 26];                                % Depth indices corresponding to available depths
    dpths_m = [1.0182, 10.537, 19.398, 29.886, 42.145, 51.38, 112.25];  % Depth values in meters
    pos = find(slice == idx);                                           % Find the position of the selected idx in the slice array
    dpth_m = dpths_m(pos);                                              % Selected depth

    slicedataT  = dataT(:, :, idx);                                     % Extract temperature at selected depth(s)
    slicedataO  = dataO(:, :, idx);                                     % Extract oxygen at selected depth(s)
    slicedataso = dataso(:, :, idx);                                    % Extract salinity at selected depth(s)

    slicedataT  = slicedataT(:);                                        % Flatten temperature matrix into a column vector
    slicedataO  = slicedataO(:);                                        % Flatten oxygen matrix into a column vector
    slicedataso = slicedataso(:);                                       % Flatten salinity matrix into a column vector

    slicedataT(slicedataT == 1e+20) = [];                               % Remove fill/missing values from temperature
    slicedataO(slicedataO == 1e+20) = [];                               % Remove fill/missing values from oxygen
    slicedataso(slicedataso == 1e+20) = [];                             % Remove fill/missing values from salinity

    slicedata = cat(2, slicedataT, slicedataso, slicedataO);            % Concatenate temperature, salinity, and oxygen into one data matrix


    %% ===================== 4. SOM commands=====================
    sData_cossarini = som_data_struct(slicedata(1:end,1:3), 'name', 'Cossarini Dataset', 'comp_names', {'Temp','Sal','O2_uM'}); % Structure
    clear sData_cossarini_newT sData_cossarini_newSO sData_cossarini_newO sData_cossarini_new                                   % Clear variables before usage
    sData_cossarini_newT(:,1)  = som_normalize(sData_cossarini.data(:,1),sData.comp_norm{2});                                   % Normalization temperature
    sData_cossarini_newSO(:,1) = som_normalize(sData_cossarini.data(:,2),sData.comp_norm{3});                                   % Normalization salinity
    sData_cossarini_newO(:,1)  = som_normalize(sData_cossarini.data(:,3),sData.comp_norm{4});                                   % Normalization oxygen


    %% ===================== 5. BEST-MATCHING UNITS (BMUs) =====================
    nRows = size(sData_cossarini_newT, 1);                                               % Number of rows in temperature dataset
    nan_col = nan(nRows, 1);                                                             % Create a column of NaN values to align dataset dimensions
    sData_cossarini_new = [nan_col, sData_cossarini_newT, sData_cossarini_newSO, sData_cossarini_newO];  % Combine NaN column + Temp + Salinity + Oxygen into one matrix

    nRowssliced = size(slicedata, 1);                                                    % Number of rows in sliced (depth-specific) dataset
    nan_colsliced = nan(nRowssliced, 1);                                                 % Create a column of NaN values for sliced dataset
    sliced_denormalize = cat(2, nan_colsliced, slicedata);                               % Concatenate NaN column and sliced data to match SOM input format

    bmus_surfaceGP = som_bmus(sMap, sData_cossarini_new);                                % Find Best Matching Units (BMUs) for surface dataset
    bmus_surfaceGP_denorm = som_bmus(sMap, sliced_denormalize);                          % Find BMUs for denormalized sliced dataset



    %% =====================================================
    %% ===== PLOTS ===== PLOTS ===== PLOTS ===== PLOTS =====
    %% =====================================================

    %% ===================== 6. Hit Plots =====================
    hit_surfaceGP = som_hits(sMap, sData_cossarini_new);                                   % Compute the number of data samples (hits) mapped to each SOM node

    % ====== Plot numeric hit counts ======
    figure(6);                                                                             % Open figure
    colormap(1-gray)                                                                       % Set colormap to inverted grayscale for better contrast
    som_show(sMap, 'norm', 'd', 'umat', 'all');                                            % Display only the U-Matrix (distance map) of the SOM
    pos = som_unit_coords(sMap);                                                           % Get the (x, y) coordinates of all SOM units (neurons)

    x_offset = 0.9;                                                                        % Horizontal offset
    y_offset = 0.9;                                                                        % Vertical offset
    pos(:,1) = pos(:,1) * sqrt(3)/1.73;                                                    % Horizontal scale
    pos(:,2) = pos(:,2) * 1.155;                                                           % Vertical scale

    for i = 1:length(hit_surfaceGP)                                                        % Loop through each SOM unit
        if hit_surfaceGP(i) > 0                                                            % Only label units that have one or more hits
            text(pos(i,1) + x_offset, pos(i,2) + y_offset, num2str(hit_surfaceGP(i)), ...  % Place text near the unit showing the hit count
                  'Color', 'r', 'FontSize', 6, 'FontWeight', 'bold');                      % Set text color, size, and bold font for visibility
        end
    end
    title(sprintf('SOM Hit Numbers Slice at -%.1f m', dpth_m));                            % Add a title to the plot
    filename = sprintf('MeHg_SOM_Hit_Numbers_Slice_at_-%.1f_m.png', dpth_m);               % Create filename based on title
    saveas(gcf, filename);                                                                 % Save current figure as PNG

    clear i k     % Clear temporary variables from workspace
end




%% =============================================================
%% ============= SAVING LOOP FOR MULTIPLE Z-LEVELS =============
%% =============================================================
z_level = [1, 5, 8, 11, 14, 16, 26];                                % Define indices of depth slices to process
slice   = [z_level];                                                % Copy z_level into slice variable
dpths_m = [1.0182, 10.537, 19.398, 29.886, 42.145, 51.38, 112.25];  % Depths in meters corresponding to indices


%% ===================== 1. Initialise =====================
baseFolder = 'Cossarini_Model_Output';                                               % Folder containing data files
fileNames = { ...                                                                    % List of .nc filenames
    '20200101_y-CMCC--PSAL-MFSe3r1-MED-b20220901_re-sv01.00.nc', ...                 % Salinity
    '20200101_y-CMCC--TEMP-MFSe3r1-MED-b20220901_re-sv01.00.nc', ...                 % Temperature
    '20200101_y-OGS--BIOL-MedBFM3-MED-b20221120_re-sv05.00.nc' ...                   % Oxygen
};
files = cellfun(@(f) fullfile(baseFolder, f), fileNames, 'UniformOutput', false);    % Create full paths for files


%% =================== 2. Load .nc files ===================
% Load once outside the loop to save time
ncid = netcdf.open(files{1}, 'NC_NOWRITE');                   % Open first file (salinity)
dataso = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'so'));    % Read salinity variable
netcdf.close(ncid);                                           % Close the file
ncid = netcdf.open(files{2}, 'NC_NOWRITE');                   % Open second file (temperature)
dataT = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'thetao')); % Read temperature variable
netcdf.close(ncid);
ncid = netcdf.open(files{3}, 'NC_NOWRITE');                   % Open third file (oxygen)
dataO = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'o2'));     % Read oxygen variable
netcdf.close(ncid);

ncid = netcdf.open(files{1}, 'NC_NOWRITE');                   % Re-open first file to get latitude
datalat = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'lat'));  % Read latitude
netcdf.close(ncid);
ncid = netcdf.open(files{1}, 'NC_NOWRITE');                   % Re-open first file to get longitude
datalon = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'lon'));  % Read longitude
netcdf.close(ncid);
ncid = netcdf.open(files{1}, 'NC_NOWRITE');                   % Re-open first file to get depth
varid = netcdf.inqVarID(ncid, 'depth');                       % Get variable ID for depth
datadepth = netcdf.getVar(ncid, varid);                       % Read depth variable

% Subset to region/depth range
dataT = dataT(12:end, :, 1:125);                   % Keep subset of temperature data
dataso = dataso(12:end, :, 1:125);                 % Keep subset of salinity data
datalon = datalon(12:end);                         % Subset longitude

% Create coordinate indexing
[lat_grid, lon_grid] = meshgrid(datalat,datalon);  % Create 2D grid of lat/lon
flat_lat = lat_grid(:);                            % Flatten latitude grid
flat_lon = lon_grid(:);                            % Flatten longitude grid
latlon_matrix = [flat_lat, flat_lon];              % Combine lat/lon into single matrix


%% ===================== 3. Loop over depth slices =====================
for n = 1:length(slice)                            % Loop over each depth slice

    %% ========== Data Preparation ===============
    latlon_sliced = latlon_matrix;                 % Copy lat/lon for current slice
    idx = slice(n);                                % Current depth index
    dpth_m = dpths_m(n);                           % Current depth in meters

    fprintf('\nProcessing depth index %d (%.1f m)...\n', idx, dpth_m);  % Print processing info

    slicedataT  = dataT(:, :, idx);                           % Extract temperature at current depth
    slicedataO  = dataO(:, :, idx);                           % Extract oxygen at current depth
    slicedataso = dataso(:, :, idx);                          % Extract salinity at current depth

    slicedataT  = slicedataT(:);                              % Flatten temperature
    slicedataO  = slicedataO(:);                              % Flatten oxygen
    slicedataso = slicedataso(:);                             % Flatten salinity

    latlon_sliced(slicedataT == 1e+20, :) = [];               % Remove lat/lon where temperature is missing
    slicedataT(slicedataT == 1e+20) = [];                     % Remove missing temperature values
    slicedataO(slicedataO == 1e+20) = [];                     % Remove missing oxygen values
    slicedataso(slicedataso == 1e+20) = [];                   % Remove missing salinity values

    slicedata = cat(2, slicedataT, slicedataso, slicedataO);  % Combine T, Sal, O2 into one matrix

    %% ========== Data Elaboration ===============
    sData_cossarini = som_data_struct(slicedata, 'name', 'Cossarini Dataset', 'comp_names', {'Temp','Sal','O2_uM'});  % Create SOM data struct

    sData_cossarini_newT  = som_normalize(sData_cossarini.data(:,1), sData.comp_norm{2});  % Normalize temperature
    sData_cossarini_newSO = som_normalize(sData_cossarini.data(:,2), sData.comp_norm{3});  % Normalize salinity
    sData_cossarini_newO  = som_normalize(sData_cossarini.data(:,3), sData.comp_norm{4});  % Normalize oxygen

    nRows = size(sData_cossarini_newT, 1);                                                               % Number of rows in normalized temperature
    nan_col = nan(nRows, 1);                                                                             % Create column of NaNs for alignment
    sData_cossarini_new = [nan_col, sData_cossarini_newT, sData_cossarini_newSO, sData_cossarini_newO];  % Combine NaN + normalized data
    bmus_surfaceGP = som_bmus(sMap, sData_cossarini_new);                                                % Compute BMUs (best matching units)
    hg_codebook = sMap.codebook(bmus_surfaceGP, 1);                                                      % Extract corresponding codebook entries
    hg_codebook_denorm = som_denormalize(hg_codebook, sMap.comp_norm(1));                                % Denormalize BMU values

    nRows = size(sData_cossarini.data, 1);                             % Number of rows in original data
    nan_col = nan(nRows, 1);                                           % Create NaN column for alignment
    sData_cossarini.data = [nan_col, sData_cossarini.data];            % Add NaN column to original data
    sData_cossarini.data(:,1) = hg_codebook_denorm;                    % Replace first column with denormalized BMU values
    sData_cossarini.data = cat(2,latlon_sliced,sData_cossarini.data);  % Collate lat/lon to data

    varName = sprintf('sData_cossarini_%dm', round(dpth_m));      % Generate variable name based on depth
    fprintf('Saving variable: %s\n', varName);                    % Print saving info
    eval([varName ' = sData_cossarini.data;']);                   % Assign data to variable with dynamic name
    outputFile = 'MeHg_Cossarini_SOM_MultiZlevel.mat';            % Output .mat file name
    save(outputFile, varName, '-append');                         % Save variable to .mat in append mode

    %% ============ Plot Hg as a Map ============
    final_lon = flip(sData_cossarini.data(:,1));    % flip longitude data
    final_lat = flip(sData_cossarini.data(:,2));    % flip latitude data
    final_hg  = flip(sData_cossarini.data(:,3));    % flip mercury values

    figure(7+n);                                                     % create a new figure for each depth
    set(gcf, 'Position', [100, 100, 1200, 400]);                     % set figure size
    scatter(final_lat, final_lon, 10, final_hg, 'filled');           % plot scatter of mercury values
    colorbar;                                                        % show colorbar
    caxis([0 0.4]);                                                  % fix color scale from 0.0 to 0.4
    xlabel('Longitude');                                             % label x-axis
    ylabel('Latitude');                                              % label y-axis
    xlim([-7 36]);                                                   % ticks x-axis
    ylim([30 46]);                                                   % ticks y-axis
    title(sprintf('Depth: %.1f m', dpth_m));                         % set title with depth
    filename = sprintf('MeHg_Cossarini_depth_%.1f.png', dpth_m);     % prepare filename
    saveas(gcf, filename);                                           % save figure as PNG

end

fprintf('\nAll depths processed and saved to %s\n', outputFile);   % Print completion message

clear i idx k n varid         % Clear temporary variables from workspace




%% =============================================================
%% ===================== ROSATI MODEL ==========================
%% =============================================================
baseFolder_rosati = 'Rosati_Data_Model_Output/';       % Base folder containing NetCDF files
fileNames_rosati = { ...                               % List of NetCDF file names to process
    'ave.20140616-00_00_00.DMHg.nc', ...
    'ave.20140616-00_00_00.MMHg.nc', ...
};
files_rosati = cellfun(@(f) fullfile(baseFolder_rosati, f), fileNames_rosati, 'UniformOutput', false);  % Create full paths for all files
allVarNames_rosati = cell(1, numel(files_rosati));                                                      % Initialize cell array to store variable names for each file

ncid = netcdf.open(files_rosati{1}, 'NC_NOWRITE');  % Open first NetCDF file (PSAL) in read-only mode
varid = netcdf.inqVarID(ncid, 'DMHg');              % Get variable ID for 'so' (salinity)
datadmhg = netcdf.getVar(ncid, varid);              % Read the 'so' variable data into MATLAB/Octave
netcdf.close(ncid);                                 % Close the first NetCDF file
ncid = netcdf.open(files_rosati{2}, 'NC_NOWRITE');  % Open first NetCDF file (PSAL) in read-only mode
varid = netcdf.inqVarID(ncid, 'MMHg');              % Get variable ID for 'so' (salinity)
datammhg = netcdf.getVar(ncid, varid);              % Read the 'so' variable data into MATLAB/Octave
netcdf.close(ncid);                                 % Close the first NetCDF file

ncid = netcdf.open(files_rosati{1}, 'NC_NOWRITE');            % Re-open first file to get latitude
datalat = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'lat'));  % Read latitude
netcdf.close(ncid);
ncid = netcdf.open(files_rosati{1}, 'NC_NOWRITE');            % Re-open first file to get longitude
datalon = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'lon'));  % Read longitude
netcdf.close(ncid);
ncid = netcdf.open(files_rosati{1}, 'NC_NOWRITE');            % Re-open first file to get depth
varid = netcdf.inqVarID(ncid, 'depth');                       % Get variable ID for depth
datadepth = netcdf.getVar(ncid, varid);                       % Read depth variable

datamehg = datadmhg + datammhg;

clear i ncid nvars varname varid     % Clear temporary variables from workspace

z_level = [1, 5, 8, 11, 14, 16, 26];                                % Define indices of depth slices to process
slice   = [z_level];                                                % Copy z_level into slice variable
dpths_m = [1.0182, 10.537, 19.398, 29.886, 42.145, 51.38, 112.25];  % Depths in meters corresponding to indices

% Create coordinate indexing
[lat_grid, lon_grid] = meshgrid(datalat,datalon);  % Create 2D grid of lat/lon
flat_lat = lat_grid(:);                            % Flatten latitude grid
flat_lon = lon_grid(:);                            % Flatten longitude grid
latlon_matrix = [flat_lat, flat_lon];              % Combine lat/lon into single matrix

%% ===================== 1. Loop over depth slices =====================
for n = 1:length(slice)         % Loop over each depth slice

    %% ========== Data Preparation ===============
    idx = slice(n);             % Current depth index
    dpth_m = dpths_m(n);        % Current depth in meters
    fprintf('\nProcessing depth index %d (%.1f m)...\n', idx, dpth_m);

    slicedataME = datamehg(:, :, idx);  % Extract combined MEHg at current depth
    slicedataME = slicedataME(:);       % Flatten MEHg

    latlon_sliced = latlon_matrix;      % Copy lat/lon for current slice


    mask = (slicedataME == 2e+20);      % Remove missing values (assuming missing values are 1e+20 as in previous code)
    latlon_sliced(mask, :) = [];
    slicedataME(mask) = [];

    sData_rosati = cat(2, latlon_sliced, slicedataME);

    %% ========== Data Elaboration ===============
    varName = sprintf('sData_rosati_%dm', round(dpth_m));
    fprintf('Saving variable: %s\n', varName);
    eval([varName ' = sData_rosati;']);

    outputFile = 'MeHg_Rosati_MultiZlevel.mat';
    save(outputFile, varName, '-append');   % Save to .mat file (append)

    %% ============ Plot MEHg as a Map ============
    final_lon = flip(latlon_sliced(:,1));   % flip longitude
    final_lat = flip(latlon_sliced(:,2));   % flip latitude
    final_mehg = flip(slicedataME);         % flip MEHg values

    figure(14+n);                                                 % new figure for each depth
    set(gcf, 'Position', [100, 100, 1200, 400]);                  % set figure size
    scatter(final_lat, final_lon, 10, final_mehg, 'filled');      % scatter plot
    colorbar;                                                     % show colorbar
    caxis([0 0.4]);                                               % fix color scale from 0.0 to 0.4
    xlabel('Longitude');                                          % label x-axis
    ylabel('Latitude');                                           % label y-axis
    xlim([-7 36]);                                                % ticks x-axis
    ylim([30 46]);                                                % ticks y-axis
    title(sprintf('ROSATI Model MEHg - Depth: %.1f m', dpth_m));  % set title with depth
    filename = sprintf('MeHg_Rosati_depth_%.1f.png', dpth_m);     % prepare filename
    saveas(gcf, filename);                                        % save figure as PNG

end

fprintf('\nAll ROSATI depths processed and saved to %s\n', outputFile);   % Print completion message

clear idx n mask latlon_sliced final_lat final_lon final_mehg    % Clear temporary variables from workspace




%% ============================================================
%% ====================== ERROR MEASURE =======================
%% ============================================================
clear  % clear workspace

%% ============ Load Data ============
load('MeHg_Cossarini_SOM_MultiZlevel.mat')
load('MeHg_Rosati_MultiZlevel.mat')
depths = {'1m', '11m', '19m', '30m', '42m', '51m', '112m'};   % list of model depth levels

%% ============ Custom Red-White-Blue Colorbar ============
function cmap = redwhiteblue(n)                 % define red-white-blue colormap
  if nargin < 1
      n = 256;                                  % default number of colors
  end
  half = floor(n/2);                            % split color range in half
  r = [(0:half-1)'/half; ones(half,1)];         % red increases from 0→1, stays 1
  g = [(0:half-1)'/half; (half-1:-1:0)'/half];  % green goes up then down
  b = [ones(half,1); (half-1:-1:0)'/half];      % blue starts at 1 then fades
  cmap = [r g b];                               % combine RGB into colormap
end


%% ============ Main Loop ============
for n = 1:length(depths)                        % Loop over each depth slice
  depth = depths{n};                            % current depth string (e.g. '1m')

  coss_name = ['sData_cossarini_' depth];       % sData_cossarini_1m
  ros_name  = ['sData_rosati_' depth];          % sData_rosati_1m

  coss_data = eval(coss_name);                  % get matrix from workspace
  ros_data  = eval(ros_name);                   % get matrix from workspace

  lat_c = coss_data(:,1);                       % Cossarini latitude
  lon_c = coss_data(:,2);                       % Cossarini longitude
  val_c = coss_data(:,3);                       % Cossarini value
  lat_r = ros_data(:,1);                        % Rosati latitude
  lon_r = ros_data(:,2);                        % Rosati longitude
  val_r = ros_data(:,3);                        % Rosati value

  coords_c = [lat_c, lon_c];                    % coordinates from Cossarini
  coords_r = [lat_r, lon_r];                    % coordinates from Rosati

  [common_coords, idx_c, idx_r] = intersect(coords_c, coords_r, 'rows'); % match by (lat, lon)

  diff_val = val_c(idx_c) - val_r(idx_r);       % difference of third column
  diff_data = [common_coords, diff_val];        % combined result [lat lon diff]

  diff_name = sprintf('diff_%s', depth);                % Create variable name (e.g. diff_1m)
  fprintf('Saving variable: %s\n', diff_name);          % Announce what is being saved
  eval([diff_name ' = diff_data;']);                    % Assign dynamically
  outputFile = 'MeHg_Difference_Cossarini_Rosati.mat';  % Save to .mat file (append mode)
  save(outputFile, diff_name, '-append');

  %% ============ Plot Difference ============
  figure(14 + n);                               % create a new figure per depth
  set(gcf, 'Position', [100, 100, 1200, 400]);  % set figure window size
  ax = gca;
  hold on;
  fill([-7 36 36 -7], [30 30 46 46], [0.8 0.8 0.8], 'EdgeColor', 'none');   % gray background (land)
  scatter(common_coords(:,2), common_coords(:,1), 10, diff_val, 'filled');  % scatter plot colored by diff
  colormap(redwhiteblue());                                                 % apply red-white-blue colormap
  colorbar;                                                                 % show color scale

  clim = max(abs(diff_val));    % find max abs difference
  caxis([-clim clim]);          % symmetric color axis limits

  xlabel('Longitude');                                                                       % X label
  ylabel('Latitude');                                                                        % Y label
  xlim([-7 36]);                                                                             % longitude range
  ylim([30 46]);                                                                             % latitude range
  title(sprintf('Cossarini - Rosati Difference at Depth: %s', depth), 'FontWeight', 'bold'); % title

  filename = sprintf('MeHg_Difference_Cossarini_minus_Rosati_%s.png', depth);  % output file name
  saveas(gcf, filename);                                                  % save figure to file

  fprintf('Depth %s: computed %d matching points and saved plot "%s".\n', depth, rows(diff_data), filename);  % print progress info

end

clear ax clim common_coords coords_c coords_r coss_data coss_name ...    % Clear temporary variables from workspace
  depth depths diff_data diff_name diff_val filename ...                 % Clear temporary variables from workspace
  idx_c idx_r lat_c lat_r lon_c lon_r ros_data n ros_name val_c val_r    % Clear temporary variables from workspace




%% ============================================================
%% ===================== CLOSE ALL PLOTS ======================
%% ============================================================
close all           # Close all plots
