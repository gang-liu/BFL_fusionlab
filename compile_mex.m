%ITK_DIR = '/Users/local_msabuncu/Documents/work/code/ITK/';

ITK_DIR = '/home/zhaolong/program/ITK3.20.0/ITK3.20.0-bin/';

if (~strcmp(comp_id(1:4), 'GLNX'))
    cd ./diffeomorphic_registration/
    compile_mex_example(ITK_DIR);
    cd ..
end
cd ./misc/toolbox_fast_marching/
compile_mex
cd ../..