clear
vol1=MRIread('label02.nii.gz');
vol2=MRIread('label12.nii.gz');
[x y z]=size(vol1.vol);
for i=1:x
    for j=1:y
        for k=1:z
            if vol1.vol(i,j,k)~=24
                vol1.vol(i,j,k)=0;
            end
            if vol2.vol(i,j,k)~=24
                vol2.vol(i,j,k)=0;
            end
        end
    end
end
    
MRIwrite(vol1,'02singlelabel_24.nii.gz');
MRIwrite(vol2,'12singlelabel_24.nii.gz');