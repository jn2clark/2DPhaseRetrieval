function [ new_array ] = zero_pad_ver3( input,newx,newy,newz )
%jclark

nd=ndims(input);

nn=size(input);

x=newx-nn(2);
y=newy-nn(1);

if ndims(input) == 3,z=newz-nn(3);else z=0;end

nnc=[floor(x/2),ceil(x/2),floor(y/2),ceil(y/2),floor(z/2),ceil(z/2)];

new_array = init_pad(input,nnc);
new_array = init_crop(new_array,nnc);

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

function [ cropped ] = crop_dim(array,new,dim )
%jclark

sz=size(array);

%switch dim
if dim == 1
    cropped=array(:,new(1):new(2),:);end
if dim == 2 
    cropped=array(new(1):new(2),:,:);end
if dim == 3
    cropped=array(:,:,new(1):new(2));end
%end

end

