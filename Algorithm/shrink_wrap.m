function support = shrink_wrap(pn,threshold,sig,type)
% - Jesse Clark, LCN, UCL October-November 2010
%   jesse.clark@ucl.ac.uk, jessenclark@gmail.com
%shrinkwrap 
%five types, 'gauss' and 'box', 'gauss_percent' and 'percent','gauss_fill'
%gauss does gauss smoothing and box does boxcar smoothing
%sig is the smoothing paramter, is equal to the standard deviation for
%gauss and is equal to the box length for box

try
    type;
catch
    type='gauss';
end

support=abs(pn)-abs(pn);

dims=ndims(support);

if dims == 2, sigx=[sig,sig];end  %this is for gauss_conv_fft
if dims == 3, sigx=[sig,sig,sig];end

if strcmp(type,'gauss_percent')
    
    if dims == 2, gauss=gauss_2D(7,7,sig,sig,0);end
    if dims == 3, gauss=gauss_3D(7,7,7,sig,sig,sig);end
    
    smooth=gauss_conv_fft(abs(pn),sigx ,0);
    
    smooth=smooth/max(max(max(smooth)));
    [h x]=hist(smooth(:),10000);
    cs=cumsum(h)/numel(smooth(:));
    ind=find(cs >= (1e0-threshold) );
    
    th=x(ind(1));
    ind=(smooth >= th );
    support=ind.*1e0;
    
end

if strcmp(type,'kmeans')
    
    ind=abs(pn(:));
    [IDX]=kmeans(ind,2);
    
    support=ones(size(ind));
    support(IDX == 1)=0;
    support=reshape(support,size(pn));
    
end

if strcmp(type,'percent')

    ind=(abs(pn) >= 0);
    [h x]=hist(abs(pn(ind)),10000);
    cs=cumsum(h)/numel(ind);
    ind=find(cs >= (1e0-threshold) );
    
    th=x(ind(1));
    ind=(abs(pn) >= th );
    support=ind.*1e0;
    
end

if strcmp(type,'gauss')
    
    if dims == 2, gauss=gauss_2D(7,7,sig,sig,0);end
    if dims == 3, gauss=gauss_3D(7,7,7,sig,sig,sig);end
    
    %smooth=convn(abs(pn),gauss,'same');
    
    smooth=gauss_conv_fft(abs(pn),sig ,0);
    
    smooth=smooth/max(max(max(smooth)));
    
    ind=(smooth >= threshold*max(max(max(smooth))) );
    support(ind)=1e0;
    ind=0;
    
end    
if strcmp(type,'gauss_fill')
    
    if dims == 2, gauss=gauss_2D(7,7,sig,sig,0);end
    if dims == 3, gauss=gauss_3D(7,7,7,sig,sig,sig);end
    
    %smooth=convn(abs(pn),gauss,'same');
    
    smooth=gauss_conv_fft(abs(pn),sigx ,0);
    
    smooth=smooth/max(max(max(smooth)));
    
    ind=(smooth >= threshold*max(max(max(smooth))) );
    support(ind)=1e0;
    ind=0;
    support=fill_support(support);
end    

if strcmp(type,'box')
    
    if dims == 2, gauss=padarray(zeros(sig,sig)+1e0,[7,7]);end
    if dims == 3, gauss=padarray(zeros(sig,sig,sig)+1e0,[3,3,3]);end
    
    smooth=convn(abs(pn),gauss,'same');
    smooth=smooth./max(max(max(smooth)));
    
    ind=(smooth >= threshold*max(max(max(smooth))) );
    support(ind)=1e0;
    
    

end


end