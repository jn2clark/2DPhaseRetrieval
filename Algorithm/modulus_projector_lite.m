function [pnm error params] = modulus_projector_lite(pn,data,params)
% jclark
% modulus constraint projector

% get estimate for scattered wave
psi = fftxy(pn,1);

% calculate the current error
error=calc_chi(abs(psi(data ~= 0)) ,sqrt(data(data ~= 0)));

% replace modulus
psi=replace_modulus(psi,sqrt(data),params);

% propagate back to sample plane
pnm=fftxy(psi,-1);        

end

function G = replace_modulus(F,data,params)
% jclark
% replace modulus of array while retaining phase

% get modulus
Mk=abs(F);

% init ratio array
ratio=Mk-Mk+1;

% check for zeros before dividing
ind=(Mk ~= 0);
ratio(ind)=(data(ind))./(Mk(ind));

% let zeros float if specified
if isempty(params.mask) ~= 1
    
    % only let values within mask region and where data is zero
    % float
    ratio(params.mask ~= 0 & data == 0) =1;
    
elseif params.no_zero == 1
    
    % if no mask specified, let all zero values float
    ratio(data == 0) = 1;
    
end

% multiply original function by the ratio of magnitudes
G=(ratio).*F;   

end

function error = calc_chi(Mk,Mm)
% jclark
% calculates error between two arrays based on intensity.

nume=sum(sum(sum(abs(Mk-Mm).^2)));    
denom=sum(sum(sum(abs(Mm).^2)));
error=nume/denom;

end
