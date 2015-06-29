function [params] = iterative_phasing_lite(params)
% jclark 
% performs phasing from diffraction.  This version
% specifially designed for 2D XFEL data with or without missing
% data.  params is a structure created using Matlab_phasing_ver1_1.m file

% set defaults if they don't exist
params = set_defaults(params);

% init the rand numb generator for reproducibility
rng(params.GA_random_seed);

% need to populate the initial guess 
% and generate a support or load previous reconstruction   
params = init_phasing_lite(params);
    
% get power for error calc
denom = sum(params.data(:));

% now iterate over generations/population/iteration. 
for pp=1:params.generations

    % do the low-to-high resolution phasing here
    if params.GA_lres == 1
        
        disp(' ')
        disp('Performing low-to-high resolution phasing....')
        disp(' ')
        
        % change the sigma of the data mask
        params.sigma = params.GA_sig_custom(pp);
        
        if pp == 1
            dataO = params.data;  %keep oringal data
        end
        
        % create the mask for low to high res phasing
        gmask=circshift(make_data_mask(params,pp,dataO),reverse(params.GA_mask_offset));
        
        % multiply with gauss mask to give lower resolution
        params.data = dataO.*gmask;
        
    end
    
    % iterate over population members. could use parfor to
    % parralize here.
    for ww=1:params.population
        
        % flip vars are for changing algorithm and turning on
        % and off the support update
        params.flip=1;
        params.shrink_flip=-1;
        params.phase_flip=-1;
        pnm_avg=0;
        
        % init the local object estimates/iterates/support
        pn =  params.pn(:,:,ww);
        pnm = params.pnm(:,:,ww);
        support=params.support(:,:,ww);
        
        % counter for averaging
        avg_count=0;
        
        for qq=1:params.iterations
        
            % keep current it no
            params.itno = qq;
            
            % check alg state, will flip sign when met
            params.flip = calc_state(params.trigger,qq).*params.flip;    

            % check if we should update support
            params.shrink_flip = calc_state(params.shrink_trigger,qq).*params.shrink_flip; %check sw
            
            % check if we are constraining the phase
            params.phase_flip = calc_state(params.phase_trigger,qq).*params.phase_flip;
            
            % check which algorithm
            params = check_ALG( params );
            
            % generate the next input
            pn = next_iterate_lite(pnm,support,pn,params);

            % do modulus projection and calc error 
            [pnm chi(qq) params] = modulus_projector_lite(pn,params.data,params);

            
            % output for diagnostics, check it no, population, generation
            % and flip state for support update
            if params.silent ~= 1
                disp([num2str(qq),', ',num2str(pp),', ',num2str(ww),', ',params.ALG,', ',num2str(chi(qq)),', ',num2str(params.shrink_flip),', ',num2str(params.phase_flip)])
            end
            
            % update support 
            if params.shrink_flip == 1 && params.threshold > 0
                
                switch params.sw_type
                    
                    case 'gauss-minarea'
                        support=shrink_wrap(abs(pnm),params.threshold,params.sigma,'gauss');
                        
                        if sum(support(:)) < params.sw_min_area*numel(support(:))
                            support=shrink_wrap(abs(pnm),params.sw_min_area,params.sigma,'gauss_percent');
                            disp(['Minimum support area not met, setting to minimum of ',num2str(params.sw_min_area)])
                        end
                    case 'gauss-maxarea'
                        support=shrink_wrap(abs(pnm),params.threshold,params.sigma,'gauss');
                        
                        if sum(support(:)) > params.sw_min_area*numel(support(:))
                            support=shrink_wrap(abs(pnm),params.sw_max_area,params.sigma,'gauss_percent');
                            disp(['Maximum support area exceeded, setting to maximum of ',num2str(params.sw_min_area)])
                        end
                        
                    case 'gauss-comp'    
                        
                        support=shrink_wrap(abs(pnm).^.5,params.threshold,params.sigma,'gauss');
                
                    otherwise
                        
                        support=shrink_wrap(pnm,params.threshold,params.sigma,params.sw_type);
                end
                
            end
            
            % do some averaging towards the end
            
            if qq > params.iterations-params.iterate_avg_its && params.iterate_avg_save == 1
               disp('Averaging....')
               pnm_avg=pnm+pnm_avg;
               avg_count = avg_count+1;
            else
               pnm_avg = pnm;
               avg_count=1;
            end
            
            
        end %end iterations loop
        
        % keep recently iterated and update structure
        params.chi(pp,ww,:)= chi;   
        params.pnm(:,:,ww) = pnm_avg/avg_count;
        params.pn(:,:,ww) = pn;
        params.support(:,:,ww) =support;
        
    end %end population loop
    
    % align the iterates and zero phase
    % only necessary if more than 1 random start
    if params.population > 1
        params.pnm = align_iterates_lite(params.pnm,1,0,0);
    end
    
    % now do the generation of new iterates
    % don't do after last generation as we want 
    % to return estimates
    if pp ~= params.generations
        if params.population > 1    
            [params] = breed_iterates_lite(params);
        end
    else
        % make an average and output as well
        if params.GA_avg_n > 0
            % average all or only some
            params.pnm_avg = mean(params.pnm(:,:,1:params.GA_avg_n),params.ndim_xyz+1);
        else
            if params.population > 1
                % cluster the outputs
                % based on image error and sharpness
                % can change this.
                params = get_pn_atts(params);
                grps = kmeans([params.chi_fin,params.sharp_fin],2);
                params.pnm_avg(:,:,1) = mean(params.pnm(:,:,(grps == 1)),params.ndim_xyz+1);
                params.pnm_avg(:,:,2) = mean(params.pnm(:,:,(grps == 2)),params.ndim_xyz+1);
            else
                params.pnm_avg(:,:,1) = params.pnm;
            end
        end
    end
    
    % update support since some will have 
    % shifted and be reflected after aligning and breeding
    % not necessary when only one random start
    if params.population > 1
        params = generate_new_supports(params);
    end
    
end %end generation loop

if params.GA_lres == 1
    %return the original data if it was masked
    params.data = dataO;  
end


end

function params = set_defaults(params)
% jclark
% sets some initial defaults if the
% fields do not exist.  always add a 
% default if a new paramater is added somewhere.
% allows backwards compatability

try
    params.generations;
catch
    params.generations=1;
end

try
    params.population;
catch
    params.population=1;
end

try 
    params.norm_to_data;
catch
    params.norm_to_data = 1;
end

try
    params.mask;
catch
    params.mask=[];
end

try 
    params.GA_mask_offset;
catch
    params.GA_mask_offset = [0,0];
end

try
    params.GA_avg_n;
    if params.GA_avg_n > params.population;
        params.GA_avg_n = params.population;
    end
catch
    params.GA_avg_n = params.population;
end

try
    params.save_images;
catch    
    params.save_images = 1;
end

try
    params.phase_range; 
catch
    params.phase_range = [-pi,pi];
end

try
    params.phase_trigger;
catch
    params.phase_trigger = -1;
end

try
    params.sw_min_area;
catch
    params.sw_min_area = 0;
end

try
    params.sw_max_area;
catch    
    params.sw_max_area = 1;
end

end

function flip=calc_state(trigger,qq)
% jclark
% returns 1 or -1 if iteration number matches one 
% of trigger elements
% if there is none returns 1 else returns -1

xx=fix(trigger)-fix(qq);
state=( xx == 0 );
state=sum(state);
if state == 1, flip = -1; else flip = 1;end;

end

function [ params ] = check_ALG( params )
% jclark
% switches the algorithm 

if params.flip == 1
    params.ALG=params.ALG1;
else
    params.ALG=params.ALG2;
end


end

function array = real_positive(array)
% jclark
% does a projection onto R+ for an array

array = real(array);

array(array <0)=0;



end

function support = shrink_wrap(pn,threshold,sig,type)
% jclark
% generate a support - binary mask for an image
% five types, 'gauss' and 'box', 'gauss_percent' and 'percent','gauss_fill'
% gauss does gauss smoothing and box does boxcar smoothing
% sig is the smoothing paramter, is equal to the standard deviation for
% gauss and is equal to the box length for box

% set a default
try
    type;
catch
    type='gauss';
end

% allocate support array
support=zeros(size(pn));

% get dimsnionality
dims=ndims(support);

% this is for gauss_conv_fft
if dims == 2, sigx=[sig,sig];end  
if dims == 3, sigx=[sig,sig,sig];end

% make gauss mask for init smoothing
if dims == 2, gauss=gauss_2D(7,7,sig,sig,0);end
if dims == 3, gauss=gauss_3D(7,7,7,sig,sig,sig);end

switch type
    
    case 'gauss_percent'
        % will set the support area to a fixed value.  Determined after smoothing
        % with a gauss kernal
        
        % save time by making gauss directly in Fourier domain
        smooth=gauss_conv_fft(abs(pn),sigx ,0);

        % normalize
        smooth=smooth/max(max(max(smooth)));

        % make histogram
        [h x]=hist(smooth(:),10000);

        % do cumuluative to find threshold
        % for a fixed percent of pixels
        cs=cumsum(h)/numel(smooth(:));
        ind=find(cs >= (1e0-threshold) );
        
        %get the threshold
        th=x(ind(1));

        % now use the threshold to get the support
        ind=(smooth >= th );
        support=ind.*1e0;
        
    case 'kmeans'
        % does a kmeans clustering to get the support
        
        % make into 1D
        ind=abs(pn(:));
        
        % cluster
        [IDX]=kmeans(ind,2);
        
        % make the support and fill
        support=ones(size(ind));
        support(IDX == 1)=0;
        support=reshape(support,size(pn));
        
    case 'percent'
        % calculates the support without prior smoothing
        % and finds the threshold that gives a fixed percent
        % of the total array
        ind=(abs(pn) >= 0);
        [h x]=hist(abs(pn(ind)),10000);
        cs=cumsum(h)/numel(ind);
        ind=find(cs >= (1e0-threshold) );

        th=x(ind(1));
        ind=(abs(pn) >= th );
        support=ind.*1e0;   
        
    case 'gauss'
        
        % calculates the support by smoothing amplitude
        % with a gauss kernal then thresholding
        
        % save time by calculating FFT of gauss analytically
        smooth=gauss_conv_fft(abs(pn),sig ,0);

        smooth=smooth/max(max(max(smooth)));

        ind=(smooth >= threshold*max(max(max(smooth))) );

        if threshold*max(max(max(smooth))) > min(smooth(:))
            disp('Support threshold too low...')
            ind=(smooth >= 2*threshold*max(max(max(smooth))) );
        end

        support(ind)=1e0;
        ind=0;
    
    case 'gauss_fill'
        
        % same as gauss but will also
        % shrink support and take intersection of all
        % yielding a filled in support
    
        smooth=gauss_conv_fft(abs(pn),sigx ,0);

        smooth=smooth/max(max(max(smooth)));

        ind=(smooth >= threshold*max(max(max(smooth))) );
        support(ind)=1e0;
        ind=0;
        support=fill_support(support);
        
end


end

function gmask = make_data_mask(params,ww,data)
% jclark
% makes a gaussian mask for the data

% get the sigma for the smoothing
sigxy = calc_scale_fact(params,ww);

% check if custom values provided, otherwise use a linear scaling
if isempty(params.GA_lres_custom) ~= 1,sigxy=params.GA_lres_custom(ww);disp('Using custom lres values....');end

% if the sigmas are greater than array size
% just return unity.  use this as the cutoff
if sigxy < 1,
    disp(['USING GAUSS MASK - ',num2str(sigxy)])
    if ndims(data) == 2,
        
        gx=size(data,2);
        gy=size(data,1);
        gmask=gauss_2D(gx,gy,gx*sigxy,gy*sigxy);
       
    end

    if ndims(data) == 3,
        gx=size(data,2);
        gy=size(data,1);
        gz=size(data,3);
        gmask=gauss_3D(gx,gy,gz,gx*sigxy,gy*sigxy,gz*sigxy);
    end
else
    gmask=1; 
end


end


function sigxy = calc_scale_fact(params,ww)
% jclark
% makes a scaling based on iterations and current iteration

sigxy=( (ww-1)/(params.GA_lres_genstop-1)*(1-params.GA_lres_init)+params.GA_lres_init).^params.GA_lres_pow;
                
if sigxy < params.GA_lres_init,sigxy=params.GA_lres_init;end

end

function pn = next_iterate_lite(pnm,support,pn,params)
% jclark
% create the next iterate

switch params.ALG

    case 'OSS'
        params.OSS_sigx(qq) = qq/params.iterations * params.nn(2)/2;
        params.OSS_sigy(qq) = qq/params.iterations * params.nn(1)/2;
        pn = pnm.*support + convn((1-support).*(pn - params.beta * pnm),gauss_2D(7,7,params.OSS_sigx(qq),params.OSS_sigy(qq)),'same');

    case 'HIO'                
        pn = phase_constraint_lite(pnm,support,params.phase_range,params.phase_flip) + (pn - params.beta * pnm) - (phase_constraint_lite(pn,support,params.phase_range,params.phase_flip) - params.beta.* phase_constraint_lite(pnm,support,params.phase_range,params.phase_flip));
    case 'HIOr'
        pn = real(pnm.*support) + (pn - params.beta * pnm) - real(support.*(pn - params.beta * pnm));
    case 'HIOrp'
        pn = real_positive(pnm.*support) + (pn - params.beta * pnm) - real_positive(support.*(pn - params.beta * pnm));
    case 'ER'
        pn = phase_constraint_lite(pnm,support,params.phase_range,params.phase_flip);
    case 'ERr'
        pn = real(pnm.*support) ;
    case 'ERrp'
        pn = real_positive(pnm.*support) ;
    case {'SF','CF'}
        pn = pnm.*support - (1-support).*(pnm);
    case {'RAAR'}
        pn = 0.5*params.beta*( (2*support-1).*(2*pnm-pn)+pn)+(1-params.beta).*pnm;
        
end
            
end