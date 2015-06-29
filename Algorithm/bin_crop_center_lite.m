function [ params ] = bin_crop_center_lite(params)
%jclark
%loads data, aligns data, centers data, crops data, bins data
%thresholds data, removes aliens (spurious data)
%returns the intensity

%set defualts
params = set_params_defaults(params);

%background subtract flag, 1 =yes.  will turn off if no file found
do_bg = 1;

%create full paths to files and check bg files
params.full_files=strcat(params.data_dir,params.files);
if numel(params.back) ~= 0,params.full_bg=strcat(params.data_dir,params.back);else params.full_bg=[];do_bg=0;end

params.nfiles = size(params.files,2);

%identify files type
tmp_str = lower(params.files{1});
params.file_type = tmp_str(end-2:end);

%get the load file function aliased
switch params.file_type
    
    case {'tif','iff'}
        fileload = @load3Dtif;
    case {'mat'}
        %fileload = @load3Dmat;
        fileload = @load3Dmat_S;
    case {'hdf'}
        fileload = @load3Dhdf;
    otherwise
        error('**Error loading data.**  Unsupported file type.  Add to bin_crop_center_lite.m')
        
end

disp('#################################')
disp('Loadong data....')

%now load the files
for qq = 1:params.nfiles
    
    %load first and get dims
    temp = fileload(params.full_files{qq});
    
    %check for bg
    %assumes bg file provided for all, can change this if necessary
    if do_bg == 1
        bg = fileload(params.full_bg{qq});
    else
        bg = 0;
    end
    
    if qq == 1    
        params.ndim_xyz = ndims(temp);
    end
    
    switch ndims(params.ndim_xyz)
        
        %will have 2 spatial dims and one other
        case 2
            data(:,:,qq) = temp - bg;
            temp=[];
            
        %will have 3 spatial dims plus one
        case 3
            data(:,:,:,qq) = temp - bg;
            temp=[];
    end
        
end

bg =[];

%disp('#################################')
disp('Done....')
disp(' ')

%now align the data and sum if we have multiple sets
%could do it in the other loop but neater outside, no penalty in
%performance for the data we use
if params.nfiles > 1
    disp('Aligning data....')
    disp(' ')
    data = sum(align_data_lite(data),ndims(data));
    disp('Done....')
    disp(' ')
end

%remove aliens (spurious scatter if it exists)
data = remove_aliens(data,params);

%now do the padding and cropping, nnc is a 6 elements vector with the 
%amount to take off each dimension, + pads, -ve crops 
data = init_pad(data,params.nnc);
data = init_crop(data,params.nnc);

%cehck fft pad size, want powers of 2 or small primes x powers of 2
if params.no_fft_pad == 0
    data=fft_pad(data,[params.binning,1]);
end

%now threshold the data
if params.subtract_dc == 0
    disp(['Values below -[',num2str(params.min_data),']- will be set to 0...'])    
    data(data < params.min_data) = 0;
else
    disp(['Values below -[',num2str(params.min_data),']- will be subtracted and then zeroed...'])
    data = data - params.min_data;
    data(data <0 ) = 0;
end
disp(' ')

%check the first photon peak
if params.no_hist == 0
    plot_hist_data(data,0)
end

%bin the data.  will pad if not a binning is not a factor of current size
[data] = bin_data_lite(data,params);

%can also do a secondary threshold after the binning
if isempty(params.schot_th) ~= 1
    disp(['Doing secondary thresholding - [',num2str(params.schot_th),']'])
    ind=(data < params.schot_th);
    data(ind)=0;
end

params.data = data;
data=[];

%get some other things from the data
params.nn=size(params.data);

if numel(params.nn) == 2, params.nn(3)=1;end


end

function params = set_params_defaults(params)

try
    params.binning;
catch
    params.binning=[1,1];
end


try
    params.no_center;
catch
    params.no_center=0; %default to centering (0=yes, 1=no)
end

try
    params.no_fft_pad;
catch
    params.no_fft_pad=0; %default to fft padding (0=yes, 1=no)
end

try
    params.pad_ptych;
catch
    params.pad_ptych=0; %default is to do 3d padding rather than 2d for ptych
end
try 
    params.center_ptych;
catch
    params.center_ptych=0;
end
try
    params.no_hist;
catch
    params.no_hist=0;
end

try
    params.bg_mult;
catch
    params.bg_mult=1;
end
try
    params.do_2D;
catch
    params.do_2D=0;
end
try
    params.subtract_dc;
catch
    params.subtract_dc=0;
end

try
    params.det_orient;
catch
    params.det_orient=[];
end
try
    params.return_orig_size;
catch
    params.return_orig_size=0;
end
try
    params.data_shift
catch
    params.data_shift=0;
end

try
    params.skipping;
catch
    params.skipping=0;
end

try
    params.remove_noise_median;
catch
    params.remove_noise_median=0;
end

try    
    params.bin_median;
catch
    params.bin_median=0;
end

try
    params.aliens_type;
catch
    params.aliens_type='box';
end

try
    params.rem_hot_pixels;
catch
    params.rem_hot_pixels=0;
end

try
    params.laser_on;
catch
    params.laser_on = 1;
end

try
    params.schot_th;
catch
    params.schot_th=0;
end

try
    params.nnc;
catch
    params.nnc=0;
end

end

function [ data] = load3Dtif( fname )
%loads a 3d tiff into array 

temp = imfinfo(fname, 'tif');
            
num_images=numel(temp);
for qq = 1:num_images

    ff0 = imread(fname,qq);
    
    %pre-allocate
    if qq == 1, data=zeros([size(ff0),num_images]);end

    data(:,:,qq) = ff0;

end



end

function [ data] = load3Dmat_S( fname )

%load the data as data or array
load(fname)

try
    data;
catch
    data=frame;
    array=[];
end


end

function [ data] = load3Dhdf( fname )

%save the data as data or array
data=hdf5read(fnames,'data');

end

function data = remove_aliens(data,params)

%leave this out unless we need it, can also provide a mask
%to remove these
%if isfield(params,'rem_hot_pixels')
%    if params.rem_hot_pixels == 1
%        data=remove_hot_pixels(data);
%    end
%end

switch lower(params.aliens_type)

    case {'box'}
        data = remove_aliens_box(data,params.aliens);
    case {'cylinder'}
        data = remove_aliens_cylinder(data,params.aliens);    
    otherwise
        data = remove_aliens_box(data,params.aliens);
    
end
        
end

function data = remove_aliens_box(data,aliens)
%% determine aliens
a_c=0;      %alien counter for the loop
disp(' ')
disp('<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>')
disp('determining alien (<>) removal....')

nx=size(data);

if sum(sum(aliens)) == 0,
    disp('no alien removal....')
else
    disp('aliens detected, preparing to remove....')
    tri=numel(aliens)/3.0;
    nal=numel(aliens)/6.0;          %6 comes from needing 2 pairs of 3 coords
    
    if mod(nal,1) == 0,             %do removal if there are the correct number of paramters
        xyz=reshape(aliens,3,nal*2);        %puts it into the format  [x0,x1,x2...]
                                            %                         [y0,y1,y2...]
                                            %                         [z0,z1,z2...]    
        for hh = 1:2:tri,               %do every second to get the pairs of points
            a_c=a_c+1;
            disp(['removing alien - [',num2str(a_c),'/',num2str(nal),']'])
            
            xx=sort(xyz(1,hh:hh+1));  %x points, automatcially gets the corect order
            yy=sort(xyz(2,hh:hh+1));  %y points  low to high
            zz=sort(xyz(3,hh:hh+1));  %z points
            
            if min(xx) == -1,xx=[max(xx),nx(2)];end
            if min(yy) == -1,yy=[max(yy),nx(1)];end
            if min(zz) == -1,zz=[max(zz),nx(3)];end
            
            
            disp(['[xstart xfinish] = [',num2str(xx),']'])
            disp(['[ystart yfinish] = [',num2str(yy),']'])
            disp(['[zstart zfinish] = [',num2str(zz),']'])
            
            data(yy(1):yy(2),xx(1):xx(2),zz(1):zz(2))=0;
            
        end                                    
    else
        disp(' ')
        disp('#<># ERROR #<># ERROR #<># ERROR #<># ERROR #<># ERROR #<>#')
        disp(' ')
        disp('incorrect number of [x,y,z] pairs.  require 6 points per alien.')
        disp('no alien removal today....')
        disp(' ')
        disp('#<># ERROR #<># ERROR #<># ERROR #<># ERROR #<># ERROR #<>#')
    end
end
disp(' ')
disp('<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>') 

end

function data = remove_aliens_cylinder(data,aliens)
%% determine aliens
%aliens should be [xcent,ycent,rad,zstart,zstop];


a_c=0;      %alien counter for the loop
disp(' ')
disp('<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>')
disp('determining alien (<>) removal....')

nx=size(data);

if sum(sum(aliens)) == 0,
    disp('no alien removal....')
else
    disp('aliens detected, preparing to remove....')
    tri=numel(aliens)/5.0; %x,y,z, rad height
    
    
    if mod(tri,1) == 0,             %do removal if there are the correct number of paramters
        xyz=reshape(aliens,5,tri);        %puts it into the format  [x0,x1,x2...]
                                            %                         [y0,y1,y2...]
                                            %                         [z0,z1,z2...]    
        for hh = 1:tri,               %do every one 
            a_c=a_c+1;
            disp(['removing alien - [',num2str(a_c),'/',num2str(tri),']'])
            
            xx=(xyz(1,hh));  %x points, automatcially gets the corect order
            yy=(xyz(2,hh));  %y points  low to high
            rr=(xyz(3,hh));  %z points
            
            zb=(xyz(4,hh));
            zt=(xyz(5,hh));
            
            circ=zero_pad_ver3(generate_circle(ceil(2*rr),rr),nx(2),nx(1)); %just make a 2d circle
            
            xs=xx-round(nx(2)/2);
            ys=yy-round(nx(1)/2);

            for gg=zb:zt
                
                circs=1-circshift(circ,[ys xs]);
                data(:,:,gg)=data(:,:,gg).*circs;
                
            end
            
            
            disp(['[xcenter,ycenter] = [',num2str(xx),',',num2str(yy),']'])
            disp(['[rad] = [',num2str(rr),']'])
            disp(['[zstart zfinish] = [',num2str(zb),',',num2str(zt),']'])
            
            
            
        end                                    
    else
        disp(' ')
        disp('#<># ERROR #<># ERROR #<># ERROR #<># ERROR #<># ERROR #<>#')
        disp(' ')
        disp('incorrect number of [x,y,z] pairs.  require 6 points per alien.')
        disp('no alien removal today....')
        disp(' ')
        disp('#<># ERROR #<># ERROR #<># ERROR #<># ERROR #<># ERROR #<>#')
    end
end
disp(' ')
disp('<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>') 

end

function data = fft_pad(data,bin)

%% set an array of values that it will automtically pad to to give 
%a good number for the FFT, i.e powers of 2 or PX2^N, where P is a low
%numbered prime i.e <7 or 9
%pref=[32,48,64,80,96,128,144,160,192,256,320,512];
pow_2=2.^(1:10);
pref=sort([2.^(4:16),2.^(3:10)*3,2.^(3:10)*5,2.^(3:10)*9]);%sort([2.^(6:10),pow_2*3,pow_2*5,pow_2*9]);

dsz=size(data);

disp('')

try
    bin;
catch
    bin=[1,1,1];
end

if numel(dsz) == 3
    sx=dsz(2);
    sy=dsz(1);
    sz=dsz(3);

    disp(['Before padding for FFT [x,y,z] - [',num2str([sx,sy,sz]),']'])
    
    while sum(sx == bin(1)*pref) == 0,sx=sx+1;end
    while sum(sy == bin(2)*pref) == 0,sy=sy+1;end
    while sum(sz == bin(3)*pref) == 0,sz=sz+1;end

    nn=[sx,sy,sz];
    disp(['Padding for FFT [x,y,z] - [',num2str([sx,sy,sz]),']'])
    
    data=zero_pad_ver3(data,nn(1),nn(2),nn(3));
end
if numel(dsz) == 2
    sx=dsz(2);
    sy=dsz(1);

    disp(['Before padding for FFT [x,y] - [',num2str([sx,sy]),']'])
    
    while sum(sx == bin(1)*pref) == 0,sx=sx+1;end
    while sum(sy == bin(2)*pref) == 0,sy=sy+1;end

    nn=[sx,sy];
    disp(['Padding for FFT [x,y] - [',num2str([sx,sy]),']'])

    data=zero_pad_ver3(data,nn(1),nn(2));
end

disp('Complete.....')
disp(' ')

%disp('----------------------------------------------------------------')


end

function plot_hist_data(data,bef_aft)


if bef_aft == 1,
    disp(' ')
    figure(90)
    disp('Plotting histogram after binning....'),
    ind=(data <10000);
    [h x]=hist(data(ind),100);
    %p1=plot(x,h);%
    p1=plot(x,log(h));
    set(p1,'Color','blue','LineWidth',1.0);
    xlabel('ADU''s');
    ylabel('frequency');
    title('data histogram');
else
    disp(' ')
    figure(89)
    disp('Plotting histogram before binning....')
    ind=(data <2000);
    [h x]=hist(data(ind),100);
    %p1=plot(x,h);%
    p1=plot(x,log(h));
    ymin=log(h(2:end));
    ymin=min(ymin(ymin>0));
    %axis([0 2000 ymin max(log(h(2:end)))])
    set(p1,'Color','blue','LineWidth',1.0);
    xlabel('ADU''s');
    ylabel('frequency');
    title('data histogram');
end

end

function data = init_crop(data,nnc) 
%jclark

%% initial cropping

if numel(nnc) == 6
    
    
    %if sum(nnc(1:2)) > 0, data=crop_dim(data,nnc(1:2),1);end
    %if sum(nnc(3:4)) > 0, data=crop_dim(data,nnc(3:4),2);end
    %if sum(nnc(5:6)) > 0, data=crop_dim(data,nnc(5:6),3);end
    
    sz=size(data);
    xs=sz(2);
    ys=sz(1);

    if ndims(data) == 3,zs=sz(3);end
    
    if sum(nnc(1:2)) < 0,
        %disp('doing intial x cropping....')
        nncc=0;
        nncc=[abs(nnc(1))+1,xs-abs(nnc(2))];
       % disp(['xs [',num2str(xs),'] --> [',num2str(1+nncc(2)-nncc(1)),']'] )
        data=crop_dim(data,nncc,1);end
    
    if sum(nnc(3:4)) < 0, 
        %disp('doing intial y cropping....')
        nncc=0;
        nncc=[abs(nnc(3))+1,ys-abs(nnc(4))];
     %   disp(['ys [',num2str(ys),'] --> [',num2str(1+nncc(2)-nncc(1)),']'] )
        data=crop_dim(data,nncc,2);end
   
    
    if sum(nnc(5:6)) < 0, 
        %disp('doing intial z cropping....')
        nncc=0;
        nncc=[abs(nnc(5))+1,zs-abs(nnc(6))];
      %  disp(['zs [',num2str(zs),'] --> [',num2str(1+nncc(2)-nncc(1)),']'] )
        data=crop_dim(data,nncc,3);end
else
    disp('if initial cropping is required')
    disp('set nnc=[-x0,-x1,-y0,-y1,-z0,-z1], where nnc gives pixels to crop')
    disp('off each dimension otherwise set nnc=[0]')
end

end

function data = init_pad(data,nnc)
%jclark

sz=size(data);
xs=sz(2);
ys=sz(1);

try
 zs=sz(3);
end
 
if numel(nnc) == 6
    
    if sum(nnc(1:2)) > 0, 
        %disp('doing intial x padding....')
        data=padarray(data,[0 abs(nnc(1)) 0],0,'pre');
        data=padarray(data,[0 abs(nnc(2)) 0],0,'post');
        %disp(['xs [',num2str(xs),'] --> [',num2str(xs+nnc(2)+nnc(1)),']'] )
    end
    
    if sum(nnc(3:4)) > 0, 
        %disp('doing intial y padding....')
        data=padarray(data,[abs(nnc(3)) 0 0],0,'pre');
        data=padarray(data,[abs(nnc(4)) 0 0],0,'post');
        %disp(['ys [',num2str(ys),'] --> [',num2str(ys+nnc(4)+nnc(3)),']'] )
    end
    
    if sum(nnc(5:6)) > 0, 
        %disp('doing intial z padding....')
        data=padarray(data,[0 0 abs(nnc(5))],0,'pre');
        data=padarray(data,[0 0 abs(nnc(6))],0,'post');
        %disp(['zs [',num2str(zs),'] --> [',num2str(zs+nnc(5)+nnc(6)),']'] )
    end
    
else
    disp('if initial pading is required')
    disp('set nnc=[+x0,+x1,+y0,+y1,+z0,+z1], where nnc gives pixels to pad')
    disp('onto each dimension otherwise set nnc=[0]')
end


end




