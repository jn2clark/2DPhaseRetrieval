function cnj = conj_reflect(array)
%jclark
F=ifftshift(fftn(fftshift(array)));

cnj=fftshift(ifftn(ifftshift(conj(F))));


end
