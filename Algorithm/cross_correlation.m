function [ c ] = cross_correlation(a,b )
%jclark


AA=ifftshift(fftn(fftshift(conj_reflect(a))));

BB=ifftshift(fftn(fftshift(b)));

CC=AA.*BB;

c=ifftshift(ifftn(fftshift(CC)));


end

