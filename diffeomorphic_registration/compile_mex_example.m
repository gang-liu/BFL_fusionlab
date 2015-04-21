function compile_mex_example(ITK_DIR)

if nargin<1
    ITK_DIR = '/home/zhaolong/program/ITK3.20.0/ITK3.20.0-bin/';
end

mex_files_cell = {'mex_deffieldharmonicenergy.cpp', ...
    'mex_deffieldjacobiandeterminant.cpp', 'mex_deffieldjacobiandist.cpp', ...
    'mex_invcondemonsforces.cpp', 'mex_velocityfieldexp.cpp', ...
    'mex_warpimage.cpp', 'mex_warplabelimage.cpp', 'mex_weightedfwdemonsforces.cpp'};

INCLUDE_FLDRS = { ITK_DIR, ...
    [ITK_DIR '/Code/Common/'], ...
    [ITK_DIR '/Code/Algorithms/'], ...
    [ITK_DIR '/Code/IO/'], ...
    [ITK_DIR '/Code/Review/'], ...
    [ITK_DIR '/Code/BasicFilters/'], ...
    [ITK_DIR '/Code/Review/Statistics/'] ...
    [ITK_DIR '/Utilities/vxl/vcl/'] ...
    [ITK_DIR '/Utilities/vxl/v3p/netlib/'] ...
    [ITK_DIR '/Utilities/vxl/vcl/iso/'] ...
    [ITK_DIR '/Utilities/vxl/vcl/gcc/'] ...
    [ITK_DIR '/Utilities/vxl/vcl/gcc-295/'] ...
    [ITK_DIR '/Utilities/vxl/vcl/gcc-libstdcxx-v3/'] ...
    [ITK_DIR '/Utilities/vxl/vcl/generic/'] ...
    [ITK_DIR '/Utilities/vxl/core/'] ...
    [ITK_DIR '/Utilities/vxl/core/vnl/'] ...
    };

%%%% To figure out the following links may take some time.... you will need
%%%% to find the appropriate extension for your system (eg, .a or .dylib)
%%%% and the location of where these files are
comp_id = computer;

if (strcmp(comp_id(1:4), 'MACI'))
    lib_ext = '3.20.0.dylib';
    LINK_DIR = [ITK_DIR '/bin/'];
    libstd_fname = '/usr/lib/libstdc++.6.0.9.dylib';
    %You probably will have to change the above two lines. The easiest way is to locate the files with names listed in LINK_LIBS below
    %LINK_DIR should point to where these files are. lib_ext should be equal to their extensions.
    % Alternatives that have worked for me are: LINK_DIR = [ITK_DIR
    % '/bin/Release/']; and lib_ext = 'a';
end

LINK_LIBS = {'libITKCommon', 'libITKIO', 'libITKStatistics', ...
    'libITKAlgorithms', ...
    'libITKBasicFilters'...
    'libITKDICOMParser'...
    'libITKIOReview' ...
    'libITKNumerics' ...
    'libitksys', 'libitkvnl', 'libitkvnl_inst', 'libitkvnl_algo', 'libitkvcl', ...
    'libitkv3p_netlib'...
   };



for m = 1:length(mex_files_cell)
    cmd = ['mex -v '  mex_files_cell{m}];
    
    for i = 1:length(INCLUDE_FLDRS)
        cmd = [cmd ' -I''' INCLUDE_FLDRS{i} ''''];
    end
    
    for i = 1:length(LINK_LIBS)
        
        cmd = [cmd ' ' LINK_DIR LINK_LIBS{i} '.' lib_ext];
        
    end
    
    cmd = [cmd ' ' libstd_fname];
    cmd = [cmd ' -output ' mex_files_cell{m}(5:end-4)];
    eval(cmd);
end
  