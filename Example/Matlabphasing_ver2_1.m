%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Matlabphasing_ver2.1.m, 2D XFEL optimised 
%2D phase retrieval - Jesse Clark, October 2010 - 2015
%                       jesclark@stanford.edu, jessenclark@gmail.com
clear
params.version='Matlab phasing version 2.1 - May 2015';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0. Determining the current directory
% leave this

params.name_of_this_file = 'Matlabphasing_ver2_1';
params.dir_file = which(params.name_of_this_file);    %finds the current folder for the phasing
params.this_dir = params.dir_file(1:findstr(params.dir_file,params.name_of_this_file)-1);   %keep just the directory

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Data preperation

params.data_dir = params.this_dir;           
                                % specify the data directory.
                                % at the moment it will create a save dir with a name that is 
                                % generated (see bottom of page) according to phasing params.  it will be
                                % saved in the directory that this phasing file is in.  this can be over-
                                % ridden by simply specifiying another save_dir.  
                                
params.save_dir = params.this_dir;

                                % the data files and background files to load.  if no bg files are needed
                                % put params.back = {}.  Accepts h5, mat and tif 
                                
                                
params.files={'data256.mat'};   % name of the data file
params.back={};                 % if no bg required, leave empty {}.  don't use {''}.
 
params.seq_let='Demo-1-';       % sequence letter, used in save ouput. can be any string
params.comments='testing';      % write comments here and this will be saved in the params file

params.binning=[1,1];           % binning for x and y

params.skipping=0 ;             % will bin (=0), will skip points (=1)

params.aliens=[0];              % aliens are spurious scatter in the data that needs to be
                                % removed
                                % set aliens=0 for no removal, otherwise input them as    
                                % aliens=[[x0,y0,z0,x1,y1,z1],[x2,y2,z2,x3,y3,z3]]
                                % will remove two instances of aliens given by the pairs
                                % #0 and #1 and another given by #2,#3. accepts as
                                % many as you like.  Points are the same as winview
                         
params.min_data=0;              % min data threshold, below this val is set to 0.  the 
                                % threshold is applied to each data set BEFORE
                                % addition and binning.
                                
params.schot_th=00;             % secondary threshold applied AFTER binning
params.subtract_dc=0;           % leave this as 0 (unless you know what you are doing), removes dc term
params.no_center=0;             % no centering of data -> no_center=1
params.no_fft_pad=0;            % no padding for fft -> no_fft_pad=1
params.no_hist=1;               % plot the histograms of the data (=0)
params.bg_mult=1;               % multiplication to apply to the bg file (if exp time is different)
                                
params.nnc=[0,0];
                                % initial cropping of data before binning.
                                % eg. nnc=[0,0,-10,-10,5,5] will do nothing to x,
                                % will crop 10 pixels off each end in y and will pad
                                % 5 pixels to each end in z.  set nnc=0 to do
                                % nothing.

                                % load the data 
[ params ] = bin_crop_center_lite(params);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Experimental Geometry
% 

params.det_px=50e-6;              % detector pixel size (m)
params.lam=.31;                   % wavelength (nm)
params.delta=0.00;                % delta (degrees)
params.gam=0;                     % gamma (degrees)
params.arm=1.5;                   % camera distance (m)

                                  % calc the real-space pixel size, leave
                                  % this
params.spx=params.arm*params.lam/params.nn(1)/params.binning(1)/params.det_px;   
params.spy=params.arm*params.lam/params.nn(2)/params.binning(2)/params.det_px;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Algorithm selection and control 

params.iterations= 501;              % Iterations to perform.

params.start_guess='auto';          
                                    % start guess for iterative procedure. 
                                    % 'flat' uses the support
                                    % 'random' uses a random object the
                                    % size of the support.  'random-data'
                                    % uses the data with a random phase.
                                    % 'auto' uses the autocorrelation and 
                                    % 'auto-rand' uses a randomly
                                    % downsampled auto correlation
                                    % if the full path to a previous
                                    % reconstruction is supplied it will
                                    % use that.  it will align the
                                    % reconstruction with the data as well.


params.norm_to_data = 1;            % normalize start guess to have same power as data
                                    % useful for missing values at beamstop
                                    
                                    % select one or two algorithms.  
                                    % Choices are 'ER','SF','HIO','RAAR',
                                    % Switching between algorithms is determined
                                    % by trigger. to use one algorithm set, ALG2=ALG1
params.ALG1='ER' ;                  % Algorithm 1.  
params.ALG2='HIO' ;                 % Algorithm 2.  

params.trigger=sort([(0:100:params.iterations)+10,(0:100:params.iterations)]);          
                                    % trigger is an array of iteration numbers 
                                    % that determine when to switch
                                    % algorithms.  ie. trigger=[10,100,110,150] will do
                                    % ALG1 until iteration 10 then switch to
                                    % ALG2 until iteration 100 then it will
                                    % switch back to ALG1 at 110 etc etc. until all iterations
                                    % are complete.  There is no limit to how many times it
                                    % can change.  The number of iterations is set by
                                    % iterations though
                
params.beta=.9;                     %HIO,RAAR,alg parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Support parameters, including shrinkwrap and phase constraint
% Parameters controlling the support.  To turn shrinkwrap OFF either set
% threshold=0 OR set the elements of shrink_trigger <0.

params.sx = 0.45*params.nn(1);                  % x,y,z support size (pixels)
params.sy = 0.35*params.nn(2);                  % you can set them relative to the array dimensions   
params.sz = 0.45*params.nn(3);                  % but you can put in pixel values as well
 
params.threshold = 0.1;               % Shrinkwrap threshold (0<threshold<1). 
                                       % threshold=0 will turn SW OFF.
params.sigma=3;                        % Sigma of gauss kernal for smoothing in shrinkwrap
params.sw_type='gauss';
                                       %'gauss' - uses a gaussian Kernal to smooth
                                       %'kmeans' - uses kmeans clustering
                                       %            to get support
                                       %'percent' - adjusts the threshold so the support is
                                       %            a constant number of voxels.  ie. .1 will make sure
                                       %            the support is 10% of the total voxels
                                       %'gauss_fill' - is the same as gauss but fills in any
                                       %               holes, basically makes the support convex
                                       %'gauss_percent' - same as percent
                                       %                  but smooths first
                                       %'percent-auto'- estimates the size
                                       %              from the
                                       %              autocorrelation and
                                       %              keeps this fixed for
                                       %              the support
                                       %'gauss-minarea' - same as gauss but
                                       %                will switch to 'gauss-percent' if
                                       %                the support volume falls below params.sw_min_area
                                       %'gauss-maxarea' - same as gauss but
                                       %                will switch to 'gauss-percent' if
                                       %                the support volume goes above params.sw_max_area

params.sw_min_area=.05;                % for use with 'gauss-minarea'
                                       % as a fraction of the total, i.e .02
                                       % is 2 percent
                                       
params.sw_max_area=.125;                 % for use with 'gauss-maxarea'
                                       % as a fraction of the total, i.e .02
                                       % is 2 percent
                                       
                              
params.shrink_trigger=generate_trigger(5,5,1900);      
                                       % the same as trigger but for controlling 
                                       % shrinkwrap.  make all <0 to keep shrinkwrap 
                                       % off or set threshold=0.
                                
params.phase_trigger = [10];          % trigger for turning on the phase 
                                       % constraint.  usage is the same as
                                       % shrink_trigger.  leave -ve to turn
                                       % off
params.phase_range = [-pi,pi];         % range to constrain the phase (if turned on)
                                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. Guided paramaters

params.generations = 1;                   % do this many generations  
params.population = 1;                    % with this many individuals

params.GA_metric = 'chi';                 % how to decide what is the 'best'
                                          % used for the breeding only
                                          % 'chi' - usual chi squared (min
                                          % chi)
                                          % 'sharpness' - {|p|^4 (min
                                          % sharpest)
                                            
                                        
params.breed_mode1 = 'sqrt_ab';           % primary breeding mode, 'none' will do no breeding
                                          % 'sqrt_ab' uses sqrt, 'mean' uses
                                          % the mean, 'rand' uses a random
                                          % scheme

params.GA_lres=0;                         % start with low res data (=1) or always use full res (=0)
params.GA_lres_init=.1;                   % initial sig, 0-1 (array size = 1)
                                
                                          % at which generation should the full res start being used
params.GA_lres_genstop = params.generations;
params.GA_lres_pow = 1;                   % how the sig scales with generation (1 - linear, 2 - quad etc)
params.GA_sig_max = 3;                    % max support sig value, will decrease linearly until params.sigma is reached

                                          % custom sig values, leave empty,[], for automatic mode
params.GA_sig_custom = make_vals(params.generations,.8,params.GA_sig_max);
                                        
params.GA_lres_custom = make_vals(params.generations,.9,params.GA_lres_init);     


%params.start_guess = '/Users/jesseclark/Documents/MATLAB/MatlabPhasing-Lite/Rec-Test-2x1-ERHIO-501-Np5-Ng1/Rec-Test-2x1-ERHIO-501-Np5-Ng1.mat';

params.GA_random_seed=1;                %seed for random number generator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7. Miscellaneous pramaters

%params.mask = (params.data == 0);        % can also define a mask.  these values will float
                                          % mask values take precedence of
                                          % no_zero flag

params.save_images = 1;                   % save the images as png?

params.silent = 0;                        % suppresses output to the screen         
params.no_zero = 0;                       % let zero measured intensity values float, 0=no,1=yes


params.iterate_avg_save = 1;              % save the average iterate
params.iterate_avg_its = 20;              % number of iteratations to average over

%% END OF USER INPUTS % END OF USER INPUTS % END OF USER INPUTS % END OF USER INPUTS

% do the phasing
params = iterative_phasing_lite(params);

% save the phasing and images
save_matlabphasing_lite(params);


%% END OF FILE %% END OF FILE 
