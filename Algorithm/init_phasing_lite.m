function params = init_phasing_lite(params)
%jclark
%init the arrays and support. need to add support for loading a support
%need to add for 3D

%create support
if numel(size(params.data)) == 3
    support=zero_pad_ver3(ones(round([params.sy,params.sx,params.sz])),params.nn(2),params.nn(1),params.nn(3) );     
else  support=zero_pad_ver3(ones(round([params.sx,params.sy])),params.nn(2),params.nn(1) );end

%preserve power
norm=sum(params.data(:));

%generate the guesses, vectorized where possible
switch params.start_guess
   
    case 'flat'
        pn = support;
        pn = repmat(support,[1 1 params.population]);
        disp('Using support as initial guess....')
        
    case {'random','rand'}
        pn=repmat(support,[1 1 params.population]).*random('uniform',0,1,[params.nn(1),params.nn(2),params.population]);
        disp('Using random support as initial guess....')
        
    case {'random-data','data-rand'}
        temp=repmat(params.data,[1 1 params.population]).*exp(i*2*pi*random('uniform',0,1,[params.nn(2),params.nn(1),params.population]));
        temp=fftxy(temp,-1);
        pn=repmat(support,[1 1 params.population]).*temp;
        temp=[];
        disp('Using data with random phase as intial guess....')
         
    case 'auto-sq'
        auto=fftxy(params.data,-1);
        pn=repmat(abs(auto).^.5,[1 1 params.population]);
        auto=[];
        
    case 'auto'
        auto=fftxy(params.data,-1);
        pn=abs(auto).^.5.*exp(i*angle(auto));
        pn(1:2:size(pn,1),:,:)=[];
        pn(:,1:2:size(pn,2),:)=[];
        
        pn=zero_pad_ver3(pn,size(auto,2),size(auto,1),1);
        pn = repmat(pn,[1 1 params.population]);
        auto=[];
        
     case {'auto-rand','rand-auto'}
        auto=fftxy(params.data,-1);
        tempp=abs(auto).^.5.*exp(i*angle(auto));
        
        sz=size(tempp);
        
        for qq=1:params.population
            
            temp=tempp;
            ind1=ceil(rand([1,sz(1)/2])*sz(1));
            temp(ind1,:)=[];
            ind2=ceil(rand([1,sz(2)/2])*sz(2));
            temp(:,ind2)=[];
            
            pn(:,:,qq)=zero_pad_ver3(temp,size(auto,2),size(auto,1));
        end
        
        auto=[];
        
    otherwise
    
        params_prev=load(params.start_guess);
        
        pn=params_prev.params.pnm;
        support=params_prev.params.support;
        params.prev =[];
        disp('Using previous reconstruction....')
        
        psi = fftxy(pn,1);
        
        %need to check alignment
        [psi] = align_data_lite(cat(3,sqrt(params.data),psi));
        
        pn = fftxy(psi(:,:,2:end),-1);
        
        disp('Done....')
end

%normalize the start guesses.  
if  params.norm_to_data == 1
    pn = conserve_power(pn,params.data);
end

params.pn=pn;pn=[];
params.support=support;support=[];

if sum(size(params.support)) ~= sum(size(params.pn))
    params.support=repmat(params.support,[1 1 params.population]);
end

params.pnm=params.pn;

end


function array = conserve_power(array,data)
%jclark
%re-norm based on a val, assumed data is amp^2

if ndims(array) ~= ndims(data)

    array = fftxy(array,1);
    
    switch ndims(array)
        
        case 3
            nume = repmat( sqrt(sum(sum(abs(array).^2,1),2)), [size(array(:,:,1)) 1]);
            denom = repmat(sqrt(sum(data(:))),size(array));
            array=array./nume.*denom;
            nume=[];denom=[];

        case 4
            
            nume = repmat( sqrt(sum(sum(sum(abs(array).^2,1),2),3)), [size(array(:,:,:,1)) 1]);
            denom = repmat(sqrt(sum(data(:))),size(array));
            array=array./nume.*denom;
            nume=[];denom=[];
            
    end
    
    array = fftxy(array,-1);

else
    
    array = fftxy(array,1);
    array=array./sqrt(sum(abs(array(:)).^2))*sqrt(sum(data(:)));
    array = fftxy(array,-1);
end


end


function [ array ] = fftxy( array ,drc)
%jclark

if exist('drc') ~= 1,drc=1;end

for qq=1:size(array,3)
    
    switch drc
        case 1
            array(:,:,qq)=fftshift(fftn(fftshift(array(:,:,qq))));
        case -1
            array(:,:,qq)=fftshift(ifftn(fftshift(array(:,:,qq))));
    end
    
end

end

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

