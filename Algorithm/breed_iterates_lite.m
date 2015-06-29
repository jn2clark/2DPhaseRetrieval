function [params] = breed_iterates_lite(params)
%jclark
%combine iterates for guided algorithm

params=set_breed_defaults(params);

%assume that each time is aligned (time here is pop)
ntimes=size(params.pnm,ndims(params.pnm));  %get the ntimes


%get the best
disp(' ')
disp(['Using ',params.GA_metric,' as the metric for breeding....'])
switch params.GA_metric
    
    case 'chi'
        [val,ind]=min(params.chi(end,:,end));
        

    case 'sharpness'
        %calc the least sharp
        %make them have the same power
        switch ndims(params.pnm)

            case 4
                temp=abs(params.pnm)./repmat(sum(sum(sum(abs(params.pnm).^2,1),2),3),[size(squeeze(params.pnm(:,:,:,1,1))) 1 1]);
                [val,ind]=min((squeeze((sum(sum(sum(abs(temp*1000).^4,1),2),3)))));
            case 3
                temp=abs(params.pnm)./repmat((sum(sum(abs(params.pnm).^2,1),2)),[size(squeeze(params.pnm(:,:,1,1))) 1 1]);
                [val,ind]=min((squeeze((sum(sum(abs(temp*1000).^4,1),2)))));
        end
        
    %case 'area'
        
        
end

params.ind_best = ind;

%now loop to get the new ones
for qq=1:size(params.pnm,ndims(params.pnm))
   switch ndims(params.pnm)
        
       case 4
           %params.pnm(:,:,:,:,qq)=sqrt(abs(params.pnm(:,:,:,:,qq)).*abs(params.pnm(:,:,:,:,ind))).*exp(0.5*i*angle(params.pnm(:,:,:,:,ind))).*exp(0.5*i*angle(params.pnm(:,:,:,:,qq)));
           alpha=params.pnm(:,:,:,ind);
           beta=params.pnm(:,:,:,qq);
           
       case 3
           %params.pnm(:,:,:,qq)=sqrt(abs(params.pnm(:,:,:,qq)).*abs(params.pnm(:,:,:,ind))).*exp(0.5*i*angle(params.pnm(:,:,:,ind))).*exp(0.5*i*angle(params.pnm(:,:,:,qq)));
           alpha=params.pnm(:,:,ind);
           beta=params.pnm(:,:,qq);
   end
   
   disp(' ')
   disp('++++++++++++++++++++++++++++++++++++++++++++++++')
   
   switch params.breed_mode
       
       case {'sqrt','sqrt_ab'}
           
           disp('Cobining using sqrt...')
            gamma=sqrt( abs( alpha).*abs(beta)).*exp(0.5*i*angle(alpha)).*exp(0.5*i*angle(beta));
       case {'sqrt_ab_pa'}
           disp('Combining using sqrt and phase alpha...')
            gamma=sqrt( abs( alpha).*abs(beta)).*exp(i*angle(alpha));
       case {'mean'}
            disp('Adding to best....');
            gamma = 0.5*(alpha + beta);
       case {'rand','random'}
            disp('Using a random mask....');
            msk = (random('uniform',0,1,[size(alpha)]) > 0.5);
            gamma = msk.*alpha +(1-msk).* beta;   
       case {'none','na','NA','None','NONE'}
            disp('No combining....')
            gamma=beta;
       otherwise
           disp('No combining....')
            gamma=beta;
   end
   
   switch ndims(params.pnm)
        
       case 4
           params.pnm(:,:,:,qq)=gamma;
           
       case 3
           params.pnm(:,:,qq)=gamma;
   end
   
   
   disp('++++++++++++++++++++++++++++++++++++++++++++++++')
   disp(' ')
   
end


end


function params = set_breed_defaults(params)


try
    params.breed_mode;
catch
    try 
        params.breed_mode1;
        params.breed_mode = params.breed_mode1;
    catch
        params.breed_mode='sqrt';
    end
end

try
    params.GA_metric;
catch
    params.GA_metric='sharpness';
end



end