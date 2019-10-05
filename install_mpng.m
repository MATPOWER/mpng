function install_mpng
% INSTALL_MPNG is a quick installer for adding all MPNG distribution files 
%   to the Matlab path.

%   MPNG Matpower - Natural Gas
%   Copyright (c) 2019 - v0.99alpha
%   Sergio García Marín - Universidad Nacional de Colombia - Sede Manizales
%   Wilson González Vanegas - Universidad Tecnológica de Pereira
%   Carlos E. Murillo Sánchez - Universidad Nacional de Colombia - Sede
%   Manizales
% 
%   This file is part of MPNG.
%   Covered by the 3-clause BSD License (see LICENSE file for details). 

%% Add current MPNG directory the the MATLAB path

fprintf('\n ----------- MPNG installation routine ----------- \n\n')

base_dir = pwd(); 
file_sep = filesep();
addpath(pwd);

%% Add internal MPNG folders and files to the MATLAB path
folders = {'Functions',...
           'Cases'...
           'Examples'
            };
for i = 1:length(folders)
    folder = folders{i};
    directory = [base_dir , file_sep, folder];
    fprintf('Adding to the path: %s\n', directory);
    addpath(genpath(directory));
end

fprintf('\n MPNG has been successfully installed! \n')

clear base_dir file_sep folder folders i directory 