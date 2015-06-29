function params = generate_new_supports(params)
%jclark

for qq = 1:size(params.support,ndims(params.support))
    
    %generate new supports
    params.support(:,:,qq) = shrink_wrap(params.pnm(:,:,qq),params.threshold,params.sigma,params.sw_type);;
    
    
end


end

