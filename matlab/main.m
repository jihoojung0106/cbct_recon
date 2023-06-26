close all; clear; clc;

%% Preparation
% Setting parameters
MyParameter;

% Getting a list of all the raw files in the directory
path = 'input/*.raw';
files = dir(path);

% Reading projection data from file
proj_size = [param.nu param.nv];

fid = fopen('airmap_fit.raw', 'rb');
airmap = fread(fid, prod(proj_size), 'uint16=>single');
fclose(fid);
airmap = reshape(airmap, proj_size);

norm_proj = zeros(proj_size(1), proj_size(2), param.nProj, 'single');
for i = 1:param.nProj
    file = fullfile(pwd, 'input', files(i).name);

    % For debugging
    fid = fopen(file, 'rb');
    if fid == -1
        error('Failed to open file: %s', file);
    end

    raw_proj = fread(fid, prod(proj_size), 'uint16=>single');
    fclose(fid);
    raw_proj = reshape(raw_proj, proj_size);

    % Normalization
    raw_proj = raw_proj ./ airmap;
    
    % Correcting data
    raw_proj(raw_proj == 0) = 1;

    raw_proj(raw_proj > 1) = 1;

    % Negative logarithm
    norm_proj(:,:,i) = -log(raw_proj);
    
end

%% FBP
% Filtering
if param.gpu == 1
    filtered_proj = gpuArray(filtering(norm_proj, param));
else
    filtered_proj = filtering(norm_proj, param);
end

% Backprojection
imgs = CTbackprojection(filtered_proj, param);

%% Saving reconstructed volume
for i=1:param.nz
    slice = imgs(:,:,i);  % initialize the slice variable
    for row = 1:size(imgs(:,:,i), 1)
        for col = 1:size(imgs(:,:,i), 2)
            % slice(row, col) = slice(row, col) * 60775.4 + 1385.80;
            if slice(row, col) < 0.022802
               slice(row, col) = slice(row, col) * 60775.4 + 1385.80;
            else
               slice(row, col) = slice(row, col) * 68036.46 + 1551.36;
            end
        end
    end
    filename = fullfile('output', sprintf('%04d.raw', i));
    fid = fopen(filename, 'wb');
    
    fwrite(fid, slice, 'int16');
    fclose(fid);
end