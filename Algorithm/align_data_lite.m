function [arrays] = align_data_lite(arrays)
%jclark
%aligns to the first

ndim =ndims(arrays);

for qq=2:size(arrays,ndim)
    
    
    switch ndim
        
        case 3
        
            arrays(:,:,qq) = align_arrays(arrays(:,:,1),arrays(:,:,qq));
            
        case 4
       
            arrays(:,:,:,qq) = align_arrays(arrays(:,:,1),arrays(:,:,qq));
    end
    
    
    
end


end

