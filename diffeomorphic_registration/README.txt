#call cmake with the following command: 
#(e.g. setenv MATLAB_DIR /autofs/space/lyon_006/pubsw/common/matlab/current/)
#(e.g. setenv ITK_DIR /autofs/cluster/con_003/users/mert/pubsw/InsightToolkit-3.20.0/)

cmake -DITK_DIR=$ITK_DIR/ -DMATLAB_ROOT=$MATLAB_DIR/ -DCMAKE_BUILD_TYPE=Release
make

#Some of the mex files can be directly called by the user (in Matlab) and are very useful
# These are:
# [def_x, def_y, def_z] = velocityfieldexp(log_def_x, log_def_y, log_def_z); 
%%% The above function "exponentiates" (i.e. integrates) the velocity field [ log_def_x, log_def_y, log_def_z] to compute the warp/deformation field [def_x, def_y, def_z] 
# warped_mov_im = warpimage(mov_im, def_x, def_y, def_z);
# warped_mov_im( isnan(warped_mov_im) ) = 0;