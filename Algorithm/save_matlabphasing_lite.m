function save_matlabphasing_lite(params)
% jclark
% saves the output from the phasing
% saves the params and images as well as
% copying the original script

disp(' ')
disp('Saving reconstruction....')
disp(' ')

% create save name from params
[ name ] = create_save_name_lite(params);

disp(name)
disp(' ')

% make the dir for the new rec
params.save_dir_final = [params.save_dir,'/',name,'/'];

% make the dir if not exist
if isdir(params.save_dir_final) ~= 1,
    mkdir(params.save_dir_final);
end

% create the full name
params.save_name = [params.save_dir_final,name,'.mat'];

% save params
save(params.save_name,'params');

disp('Done....')
disp(' ')
disp('Saving png images...')

% get the chi 'rank'
[val ind]=sort(params.chi(end,:,end));
rnk =0;

if params.save_images == 1
    % output pictures
    for qq = ind;%1:params.population

        rnk = rnk + 1;

        fh = figure ; % returns the handle to the figure object
        %subplot(2,1,1)
        set(fh, 'color', 'white'); % sets the color to white 
        imagesc(abs(params.pnm(:,:,qq)));
        axis equal tight
        title(name)
        %title(num2str(params.chi(end,qq,end)))

        print(fh, '-dpng','-r300', [params.save_dir_final,'Amp-',num2str(rnk)]);%,'[',num2str(qq),']']);
        %subplot(2,1,2)
        imagesc(abs(params.support(:,:,qq)));
        axis equal tight
        title(name)

        print(fh, '-dpng','-r300', [params.save_dir_final,'Sup-',num2str(rnk)]);%,'[',num2str(qq),']']);

    end

    %output average
    if isfield(params,'pnm_avg')
        fh = figure ; % returns the handle to the figure object
        %subplot(2,1,1)
        set(fh, 'color', 'white'); % sets the color to white 
        imagesc(abs(params.pnm_avg(:,:,1)));
        axis equal tight
        title(name)
        print(fh, '-dpng','-r300', [params.save_dir_final,'Amp-avg1']);
        
        if size(params.pnm_avg,3) > 1
            imagesc(abs(params.pnm_avg(:,:,2)));
            axis equal tight
            title(name)
            print(fh, '-dpng','-r300', [params.save_dir_final,'Amp-avg2']);
        end
    end

    if isfield(params,'data')
        fh = figure ; % returns the handle to the figure object
        %subplot(2,1,1)
        set(fh, 'color', 'white'); % sets the color to white 
        imagesc(abs(params.data).^.125);
        axis equal tight
        title(name)
        print(fh, '-dpng','-r300', [params.save_dir_final,'Data']);

        fh = figure ; % returns the handle to the figure object
        %subplot(2,1,1)
        set(fh, 'color', 'white'); % sets the color to white 
        imagesc(abs(fftxy(params.pnm_avg(:,:,1))).^.25);
        axis equal tight
        title(name)
        print(fh, '-dpng','-r300', [params.save_dir_final,'Data-R1']);

        if size(params.pnm_avg,3) > 1
            imagesc(abs(fftxy(params.pnm_avg(:,:,2))).^.25);
            axis equal tight
            title(name)
            print(fh, '-dpng','-r300', [params.save_dir_final,'Data-R2']);
        end
        
    end
end

% do a try here to make it backwards compatable
try
    sys_str =  ['cp ',params.this_file,' ',params.save_dir_final,params.name_of_this_file,'.ran'];
    system(sys_str);
end

disp(' ')
disp('Done...')

close all

end

function [ name ] = create_save_name_lite(params)
%jclark


fnames.one=['Rec','-'];

fnames.two=[params.seq_let,'-'];

fnames.three=[num2str(params.binning(1)),'x',num2str(params.binning(2)),'-'];%,num2str(params.ndim_xyz),'D-'];

fnames.four=[params.ALG1,params.ALG2,'-'];

fnames.five=[num2str(params.iterations),'-'];

fnames.six=['Np',num2str(params.population),'-Ng',num2str(params.generations)];

%fnames.seven=['T',num2str(params.ntime)];
switch params.GA_metric
    case 'chi'
        str = 'C';
    case 'sharpness'
        str = 'S';
    otherwise
        str = 'N';
end
fnames.seven = ['-',str];

if params.GA_lres == 1,
    str = 'lr';
else
    str = 'nr';
end
fnames.eight = ['-',str];

name=[fnames.one,fnames.two,fnames.three,fnames.four,fnames.five,fnames.six,fnames.seven,fnames.eight];%

end





















