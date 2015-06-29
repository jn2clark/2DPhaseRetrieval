function [ pnm ] = phase_constraint_lite(pnm,support,phase_range,phase_flip)
% jclark
% enforces a phase constraint (along with support)
% we won't center the phase values first

% if constraint it snot turned on, do not do anything
if phase_flip == 1
    
    % get the min and mac allowable phase
    phi_max=max(phase_range);
    phi_min=min(phase_range);

    % get phase and amp
    phase=angle(pnm);
    amp=abs(pnm);

    % project all points less than phase min onto the complex wedge of
    % values. amplitude is modulated for a correct projection
    ind=( phase <= phi_min);
    amp(ind)=0;%abs(cos( (phase(ind)-phi_min) )).*amp(ind);
    phase(ind)=phi_min;

    % now do for points greater the phase max
    ind=( phase >= phi_max);

    amp(ind)=0;%abs(cos( (phase(ind)-phi_max) )).*amp(ind);
    phase(ind)=phi_max;
    
    % get the final version
    pnm=support.*amp.*exp(i*phase);

else
    
   %just do regular support constraint
   pnm = pnm.*support;
    
end




end

