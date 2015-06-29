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

