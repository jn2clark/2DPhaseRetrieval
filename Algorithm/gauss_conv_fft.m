function [ convag ] = gauss_conv_fft(arrx,sig,angle )
%UNTITLED Summary of this function goes here
%   does a convolution with a gaussian kernal
% using fft's.  creates the gaussian in recipricoal space to 
%save doing another fft.  should be quick and usable with the gpu


tot_arr=sum(abs(arrx(:)));  %needed to preserve 

arrk=ifftshift(fftn(fftshift(arrx)));   %n dimensioanl fft

sz=size(arrx);  %get dims

if numel(sz) == 3, nn=[sz(2),sz(1),sz(3)];end  %put them into xyz order for gaussian creation

if numel(sz) == 2, nn=[sz(2),sz(1)];end

sigk=nn./2.0/pi./sig    ;       %get inverse guss sigmas

if numel(sz) == 2, gaussk=gauss_2D(nn(1),nn(2),sigk(1),sigk(2),angle);end
if numel(sz) == 3, gaussk=gauss_3D(nn(1),nn(2),nn(3),sigk(1),sigk(2),sigk(3));end



convag=ifftshift(ifftn(fftshift(arrk.*gaussk)));

convag=real(convag); %keep only real part

ind=(convag < 0);       %positive vaues
convag(ind)=0;

convag=convag*tot_arr/sum(sum(sum(convag)));

if sum(sig) == 0,convag=arrx;end

end

