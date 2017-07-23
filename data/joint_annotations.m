clear all
close all
clc

main_path = '/home/luisitobarcito/Documents/CNEL_shared/MultiShape/Code/'; 
annotations_path = [main_path, 'data/face_data/'];
assert(isdir(annotations_path), 'Invalid directory')
dir_info = dir(annotations_path);
filenames = {dir_info.name};

%% process filenames to keep only relevant files
file_pattern = 'bmp.mat'; 
annotation_files = {};
faces_data = {};
for iFl = 1:length(filenames)
    if regexp(filenames{iFl}, file_pattern)
        annotation_files = cat(1, annotation_files, filenames{iFl});
        tmp_data = load(strcat(annotations_path, filenames{iFl}), 'annotations');
        faces_data = cat(1, faces_data, tmp_data.annotations);
    end
end

