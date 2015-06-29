function [data] = bin_data_lite(data,params)

bin = params.binning;
nx=size(data);

%pad the array so that it bins exactly
x0=nx(2);
y0=nx(1);

if max(size(nx)) == 2,nx=[nx,1];end

disp(' ')
disp('Resizing data....')
disp(['Current data size [x,y,z] - [',num2str([x0,y0,nx(3)]),']'])

while mod(x0,bin(1)) ~=0,x0=x0+1;end
while mod(y0,bin(2)) ~=0,y0=y0+1;end


data=padarray(data,[(y0-nx(1)),(x0-nx(2)),0],0,'pre');

data_new=zeros([y0/bin(2),x0/bin(1),nx(3)]);

for qq = 1:nx(3),     
    if bin(1) == 1 && bin(2) == 1  %check if binning set
        data_new=data;
        if qq == 1,disp('No binning since bx = 1 and by = 1'),end
    else %if it is check if it is skipping or binning
        if params.skipping == 1
            if qq == 1,disp('SKIPPING POINTS');end
            data_new(:,:,qq)=imresize(data(:,:,qq),[size(data_new,1),size(data_new,2)],'nearest');
        else            
            if qq == 1,disp('BINNING POINTS');end
            data_new(:,:,qq)=box_interp(data(:,:,qq),bin(1),bin(2),params.bin_median);
        end
    end
end

data=data_new;
data_new=[];

disp(['Array size after binning [x,y,z] - [',num2str([x0/bin(1),y0/bin(2),nx(3)]),']'])


end

function [ barray ] = box_interp( array,bx,by,med )
%jclark

try
    med;
catch
    med=0;
end

sz=size(array);

x0=sz(2);
y0=sz(1);


while mod(x0,bx) ~=0,x0=x0+1;end
while mod(y0,by) ~=0,y0=y0+1;end


array=padarray(array,[(y0-sz(1)),(x0-sz(2)),0],0,'pre');

newx=x0/bx;
newy=y0/by;

barray=zeros([by,bx]);

x00=[1:bx:x0,x0+1];
y00=[1:by:y0,y0+1];

if bx == 1 && by == 1
    barray=array;
    %disp(['No binning since bx = ',num2str(bx),' & by = ',num2str(by)])
else
    for xx = 1:newx
        for yy = 1:newy

            if med == 1
                temp=(( array(y00(yy):y00(yy+1)-1, x00(xx):x00(xx+1)-1  )));
                barray(yy,xx)=median(temp(:));
            else
                barray(yy,xx)=sum(sum( array(y00(yy):y00(yy+1)-1, x00(xx):x00(xx+1)-1  )));
            end
        end
    end
end

end