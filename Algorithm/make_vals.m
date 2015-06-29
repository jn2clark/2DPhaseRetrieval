function [ array ] = make_vals(n,maxv,minv)
%jclark

if n > 1
    array=(0:(n-1))/(n-1);

    array=array*(maxv-minv)+minv;

else
    
   array=[maxv]; 
    
end
    
end

