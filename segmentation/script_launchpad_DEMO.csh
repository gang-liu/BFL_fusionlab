if ($?LD_LIBRARY_PATH) then
       setenv LD_LIBRARY_PATH "/usr/pubsw/common/matlab/7.1/bin/glnxa64/":"$LD_LIBRARY_PATH"
       setenv LD_LIBRARY_PATH "/usr/pubsw/common/matlab/current/bin/glnxa64/":"$LD_LIBRARY_PATH"
else
       setenv LD_LIBRARY_PATH "/usr/pubsw/common/matlab/7.1/bin/glnxa64/"
       setenv LD_LIBRARY_PATH "/usr/pubsw/common/matlab/current/bin/glnxa64/"
endif


setenv LD_LIBRARY_PATH "/usr/local/matlab/extern/lib/glnxa64/":"$LD_LIBRARY_PATH"
setenv LD_LIBRARY_PATH "/usr/pubsw/common/matlab/7.4/sys/os/glnxa64/":"$LD_LIBRARY_PATH"

set s = 1

while ($s <= 10)

   foreach sig (1 5 10 50)	
	
	pbsubmit -l nodes=1:ppn=2 -c "BFL_labelfusion_DEMO_MGH_wrapper $s $sig"
	sleep 10
   end
   @ s = $s + 1
end


