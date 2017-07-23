shape_code_dir = [pwd,filesep, 'shape'];
misc_code_dir = [pwd,filesep, 'misc'];
cdfHC_code_dir = genpath([pwd, filesep, 'external/cdfHC']); 

addpath(shape_code_dir);
fprintf('Added %s to the search path \n', shape_code_dir);

addpath(misc_code_dir);
fprintf('Added %s to the search path \n', misc_code_dir);

addpath(cdfHC_code_dir);
fprintf('Added %s to the search path \n', cdfHC_code_dir);

