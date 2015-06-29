function params = get_pn_atts(params)


[params.chi_fin]=(params.chi(end,:,end))';

switch ndims(params.pnm)

    case 4
        temp=abs(params.pnm)./repmat(sum(sum(sum(abs(params.pnm).^2,1),2),3),[size(squeeze(params.pnm(:,:,:,1,1))) 1 1]);
        [params.sharp_fin]=((squeeze((sum(sum(sum(abs(temp*1000).^4,1),2),3)))));
    case 3
        temp=abs(params.pnm)./repmat((sum(sum(abs(params.pnm).^2,1),2)),[size(squeeze(params.pnm(:,:,1,1))) 1 1]);
        [params.sharp_fin]=((squeeze((sum(sum(abs(temp*1000).^4,1),2)))));
end




end