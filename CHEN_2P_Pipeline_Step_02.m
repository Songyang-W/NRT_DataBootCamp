%% //////////////////////// CHEN_2P_Pipeline_Step_02 //////////////////////
% Purpose
%         - Calcium footprint/signal extraction by Constrained Non-negative
%           Matrix Factorization;
%           ref. https://github.com/flatironinstitute/CaImAn-MATLAB

% History
% 11.22.17 - integrated by JNi; pls. keep tracking any changes made
%            by any user
% 12.04.17 - added 'cellalgo' in 'Ca' specifing algorithms used for
%            extracting the ROIs;  e.g. 'CNMF' or 'M'
% 12.08.17 - updated restore-after-stop by JNi
% 12.10.17 - added delete .h5 file by JNi
% 12.11.17 - added check trial_info by JNi
% 12.13.17 - updated option for checking/reuse existed ds_data.mat by JNi
%            added data.Ts and data.sessionfoldername;
%            updated update_spatial_components_ni.m such that Y_ds.mat will
%            be stored in a sessionfolder by setting options in
%            CNMFSetParms.m; this avoids conflicts from different users

function CHEN_2P_Pipeline_Step_02(pathdirectory,anm, sessionNo, area_channel, ref_trial, ...
    sampling_rate, downsample_f, redo_trial_alignment, multi_plane, subsession, slackID)
% tic;
if nargin == 0
    pathdirectory =  'Z:\Projects\Sensorimotor\Animals\';
    anm           = 'sm002';                         % name of the animal
    sessionNo     = '14';                             % session no., e.g., '2'
     sessionName   = [anm , '-', sessionNo];         % e.g.'jn018-1'; full session name
    sessionfile   = [sessionName '.mat'];            % file destination for save
    area_channel  = 'A0_Ch1';                        % selecting areas and channel
    ref_trial     = '';       % reference template for cross-trial alignment 
                                                     % auto = '' 
                                                     % manual = 'Avg_A1_Ch0_12-56-47.tif'                                                     
    sampling_rate = 30.2;                       % Hz, full frame acquisition rate
    downsample_f  = 10;                              % integer, downsampling factor
    redo_trial_alignment = 1;                        % redo trial by trial alignment
    multi_plane = 0;
    
    subsession = '';

else
    sessionName = [anm , '-', sessionNo];               % e.g.'jn018-1'; full session name
    sessionfile = [sessionName '.mat'];                 % final session file name for save
end

% /// data & toolbox directory


if ~exist('sampling_rate','var'), sampling_rate = 32.5868; end
if ~exist('downsample_f','var'),  downsample_f = 10;       end


sessionfoldername = [pathdirectory anm filesep '2P' filesep sessionName filesep];
savefoldername = [pathdirectory anm filesep];

%handling folder directory for subsessions
if ~isempty(subsession)
    sessionfoldername = [pathdirectory anm filesep '2P' filesep sessionName filesep anm subsession '-' sessionNo filesep];
end
if~isempty(subsession)
    savefoldername = [pathdirectory anm filesep '2P' filesep anm '-' sessionNo filesep];
end
%% Provide the following information
if ispc
    analysis_path = 'Z:\Dropbox\Dropbox\Chen Lab Team Folder\Analysis Suite\PIPELINE\2P\Flip\';
    addpath(genpath(analysis_path));
    analysis_path = 'Z:\Dropbox\Dropbox\Chen Lab Team Folder\Analysis Suite\PIPELINE\2P\NoRMCorre-master\';
    addpath(genpath(analysis_path));
    analysis_path = 'Z:\Dropbox\Dropbox\Chen Lab Team Folder\Analysis Suite\PIPELINE\2P\CaImAn-MATLAB-master\'; % Should be CaImAn
    addpath(genpath(analysis_path));
    addpath('Z:\Dropbox\Dropbox\Chen Lab Team Folder\Analysis Suite\PIPELINE\SlackMATLAB\');
elseif isunix
    poolobj = gcp('nocreate');
    delete(poolobj);
    NUM_CPUS = str2num(getenv('NSLOTS'));
    N = maxNumCompThreads(NUM_CPUS);
    analysis_path = './Flip/';
    addpath(genpath(analysis_path));
    analysis_path = './NoRMCorre-master/';
    addpath(genpath(analysis_path));
    analysis_path = './CaImAn-MATLAB-master/';
    addpath(genpath(analysis_path));
    addpath('/net/claustrum/mnt/data/Dropbox/Chen Lab Dropbox/Chen Lab Team Folder/Analysis Suite/PIPELINE/SlackMATLAB/');
    parpool(NUM_CPUS);
end
%% load files
if ~isempty(subsession)
    sessionfile = [anm subsession '-' sessionNo '.mat'];
end


if isempty(findstr(area_channel, 'A0')) == 0
    load([savefoldername sessionfile], 'CaA0');
    Ca = CaA0;
elseif isempty(findstr(area_channel, 'A1')) == 0
    load([savefoldername sessionfile], 'CaA1');
    Ca = CaA1;
elseif isempty(findstr(area_channel, 'A2')) == 0
    load([savefoldername sessionfile], 'CaA2');
    Ca = CaA2;
elseif isempty(findstr(area_channel, 'A3')) == 0
    load([savefoldername sessionfile], 'CaA3');
    Ca = CaA3;
end

FOV = Ca.FOV;


% if redo_trial_alignment == 0 
%    ref_trial = Ca.ref_trial; 
% end


files = dir(fullfile([sessionfoldername 'PreProcess' filesep area_channel filesep 'Avg' filesep ],['*' area_channel '*.tif']));   % list of Avg files


% here I check for previous runs and if ds data is available we'll use that
% -- to skip this simply delete ds_data.mat
ds_data_flag = isfile([sessionfoldername 'PreProcess' filesep area_channel filesep 'ds_data.mat']);
if ds_data_flag
    disp('loading ds_data.mat')
    load([sessionfoldername 'PreProcess' filesep area_channel filesep 'ds_data.mat']);
end
%% auto detect reference trial
    

%if ref_trial is still empty after loading ds_data then we should run this,
%doesn't matter whether ds_data_flag is true or not
if isempty(ref_trial)% && ~ds_data_flag
    if exist('trial_info','var')
        if length(trial_info) == length(files)
            Ca.trial_info = trial_info;
        end
    else
        % /// read trial info from raw data
        
        % scan the middle 50% of trials
        %commented a bunch out (Mitch)
        ref_metrics = zeros(length(files),1);
        parfor i = round(length(files)/4):round(length(files)*3/4)
            disp(['potential ref trial: ', num2str(i)])
            filename   = files(i).name;
            %reg_image = imread([sessionfoldername 'PreProcess' filesep area_channel filesep 'Avg' filesep filename]);
            filename   = filename(5:end-4);
            foldername = [sessionfoldername 'PreProcess' filesep area_channel filesep ];
            
            %only load variable we care about, speeds things up ~100x
            s = load([foldername filename '.mat'],'motion_metric');
            %motion_metric = s.motion_metric;
            ref_metrics(i,1) = 1-nanstd(s.motion_metric);
%             ref_metrics_{i}{1} = nanmean(reg_image(:));
%             ref_metrics_{i}{2} = nanmean(s.motion_metric);
%             ref_metrics_{i}{3} = 1-nanstd(s.motion_metric);
        end
        
%         ref_metrics = zeros(length(files),3);
%         for i = round(length(files)/4):round(length(files)*3/4)
%             for j = 1:3
%                 ref_metrics(i,j) = ref_metrics_{i}{j};
%             end
%         end
        %files_metrics = files;
        ref_metrics(ref_metrics==0) = NaN;
%         files_metrics(ref_metrics(:,3)<prctile(ref_metrics(:,3),50)) = [];
%         ref_metrics(ref_metrics(:,3)<prctile(ref_metrics(:,3),50),:) = [];
%         [i s] = max(ref_metrics(:,3));
        %files_metrics(ref_metrics<prctile(ref_metrics,50)) = [];
        %ref_metrics(ref_metrics<prctile(ref_metrics,50)) = [];
        [~, i] = max(ref_metrics);
        %ref_trial = files_metrics(s).name; % select trial with least motion
        ref_trial = files(i).name; % select trial with least motion
        disp(['chosen ref trial: ', ref_trial])
    end
end
%%
if ~exist('ref_trial', 'var'), ref_trial = files(1).name;  end
reference_image = im2double(imread([sessionfoldername 'PreProcess' filesep area_channel filesep 'Avg' filesep ref_trial]));
if  ~ds_data_flag
cd([sessionfoldername 'PreProcess' filesep area_channel filesep ]);

% % ////////////////////// Check if alignment has completed: added by Jianguang
files_done = dir(fullfile([sessionfoldername 'PreProcess' filesep area_channel filesep ],[area_channel '*.tif']));
if redo_trial_alignment == 1
    folder_d = {};
elseif numel(files_done)
    fnm  = {files_done.name};
    folder_d = cell(size(fnm));
    for i = 1:length(fnm)
        folder_d{i} = ['Avg_' fnm{i}(1:end-4)];
    end
else
    folder_d = {};
end

% /////////////////////// Trial-by-Trial alignment
parObj = gcp('nocreate');                                 % check current Parellel Pool
% Note: when using Chen-lab cluster, make sure to check if Parellel Pool
% was enabled; otherwise it may wait forever in Matlab2016a ...

if ~isempty(parObj)      % do parallel version
    parfor i = 1:length(files)
        
        fullname = files(i).name;
        % /// check files in the job list
        flag_done = 0;
        if ~isempty(folder_d)
            flag_done = any(ismember(folder_d , fullname(1:end-4)));
        end
        
        if ~flag_done
            
            [fileloc, mat_file, time_stamp, motion_metric] = trial_align_batch(fullname, sessionfoldername, area_channel, reference_image, redo_trial_alignment, multi_plane);
            trial_info{i}.fileloc                          = fileloc;
            trial_info{i}.mat_file                         = mat_file;
            trial_info{i}.time_stamp                       = time_stamp;
            trial_info{i}.motion_metric                    = motion_metric;
            
        end
    end
    
else   % non-parallel version
    
    for i = 1:length(files)
        fullname  = files(i).name;
        % /// check files in the job list
        flag_done = 0;
        if ~isempty(folder_d)
            flag_done = any(ismember(folder_d , fullname(1:end-4)));
        end
        
        if ~flag_done
            
            [fileloc, mat_file, time_stamp, motion_metric] = trial_align_batch(fullname, sessionfoldername, area_channel, reference_image, redo_trial_alignment, multi_plane);
            trial_info{i}.fileloc                          = fileloc;
            trial_info{i}.mat_file                         = mat_file;
            trial_info{i}.time_stamp                       = time_stamp;
            trial_info{i}.motion_metric                    = motion_metric;
            
        end
    end
end
end

%% CNMF calcium extraction : downsampling files
if ispc
    ds_filename = [sessionfoldername 'PreProcess' filesep area_channel filesep 'ds_data.mat'];
elseif isunix
    ds_filename = ['/scratch/chen_step2_prh/' sessionName '/PreProcess' filesep area_channel filesep 'ds_data.mat'];
end

if ~ds_data_flag

files = dir(fullfile([sessionfoldername 'PreProcess' filesep area_channel filesep 'Avg' filesep ],['*' area_channel '*.tif']));   % list of Avg files
numFiles    = length(files);

fr          = sampling_rate;                         % frame rate
tsub        = downsample_f;                          % degree of downsampling (for 30Hz imaging rate you can try also larger, e.g. 8-10)
data_type   = 'uint16';
data        = matfile(ds_filename,'Writable',true);  % create an object which don't load data into memory

if ~isprop(data, 'F_dark') || ~isprop(data, 'Ts')    % only complete data can be reused
    
    % ///
    data.sessionfoldername = sessionfoldername;          % added by Jianguang
    data.Y      = zeros([FOV,0],data_type);
    data.Yr     = zeros([prod(FOV),0],data_type);
    data.sizY   = [FOV,0];
    F_dark      = Inf;                                % dark fluorescence (min of all data)
    batch_size  = 2000;                               % read chunks of that size
    batch_size  = round(batch_size/tsub)*tsub;        % make sure batch_size is divisble by tsub
    Ts          = zeros(numFiles,1);                  % store length of each file
    cnt         = 0;                                  % number of frames processed so far
    tt1         = tic;
    
    for i = 1:numFiles
        fullname          = files(i).name;            % added by Jerry
        [~,file_name,ext] = fileparts(fullname);      % added by Jerry
        file_name         = file_name(5:end);
        folder_name       = [sessionfoldername 'PreProcess' filesep area_channel filesep];
        %         convert_file([folder_name file_name ext],'h5',fullfile(folder_name,[file_name,'_mc.h5'])); % added by Jerry
        %         name              = fullfile(folder_name,[file_name,'_mc.h5']); % added by Jerry
        %         info              = h5info(name);             % added by Jerry
        %edited by Mitch, no need to create h5 files
        name              = fullfile(folder_name,[file_name,ext]);
        info              = imfinfo(name);  %added by Mitch, imfinfo gives tiff info in a struct with one field per frame
        if length(info)>1
            dims          = [info(1).Height,info(1).Width,length(info)];
        else
            dims          = [info.Height,info.Width];
        end
        ndimsY            = length(dims);             % number of dimensions (data array might be already reshaped)
        Ts(i)             = dims(end);
        Ysub              = zeros(FOV(1),FOV(2),floor(Ts(i)/tsub),data_type);
        data.Y(FOV(1),FOV(2),sum(floor(Ts/tsub))) = zeros(1,data_type);
        data.Yr(prod(FOV),sum(floor(Ts/tsub)))    = zeros(1,data_type);
        cnt_sub = 0;
        for t = 1:batch_size:Ts(i)
%             Y       = bigread2(name,t,min(batch_size,Ts(i)-t+1));

            infoTemp.filename = name; Y = uint16(UCLA_ReadTiff(infoTemp));
            F_dark  = min(nanmin(Y(:)),F_dark);
            ln      = size(Y,ndimsY);
            Y       = reshape(Y,[FOV,ln]);
            Y       = cast(downsample_data(Y,'time',tsub),data_type);
            ln      = size(Y,3);
            Ysub(:,:,cnt_sub+1:cnt_sub+ln) = Y;
            cnt_sub = cnt_sub + ln;
        end
        disp(['downsampled trial: ', num2str(i), ' out of: ', num2str(numFiles)])
        data.Y(:,:,cnt+1:cnt+cnt_sub) = Ysub;               % /// FIXME : A potential bug
        data.Yr(:,cnt+1:cnt+cnt_sub)  = reshape(Ysub,[],cnt_sub);
        
        toc(tt1);
        cnt                           = cnt + cnt_sub;
        data.sizY(1,3)                = cnt;
    end
    data.F_dark = F_dark;
    data.Ts     = Ts;                                       % added by Jianguang
    
end
sizY       = data.sizY;                 % size of data matrix
end
%% now run new CNMF on patches on the downsampled file, set parameters first 
patch_size = [40,40];%/// [32, 32]       % size of each patch along each dimension (optional, default: [32,32])
overlap = [8,8];  % ////  [6,6]             % amount of overlap in each dimension (optional, default: [4,4])

patches = construct_patches(sizY(1:end-1),patch_size,overlap);
K = 4;                  % number of components to be found /// 10
tau = 2;                 % std of gaussian kernel (size of neuron)//// 7 
p = 0;                   % order of autoregressive system //// 2 (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.8;         % merging threshold

options = CNMFSetParms(...
    'd1',sizY(1),'d2',sizY(2),...
    'nb',1,...                                  % number of background components per patch
    'gnb',3,...                                 % number of global background components
    'ssub',2,...
    'tsub',4,...
    'p',p,...                                   % order of AR dynamics
    'merge_thr',merge_thr,...                   % merging threshold
    'gSig',tau,... 
    'spatial_method','regularized',...
    'cnn_thr',0.2,...
    'patch_space_thresh',0.25,...
    'create_memmap', true,...
    'max_size_thr', 500,...
    'min_size_thr', 100,...
    'space_thresh', .5,...
    'min_SNR',2);
options.foldername = [sessionfoldername 'PreProcess' filesep area_channel filesep];
%%

if  ~ds_data_flag
    % Run on patches (around 15 minutes)
    % keyboard
    [A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(data,K,patches,tau,p,options);

    % compute correlation image on a small sample of the data (optional - for visualization purposes)
    Cn = correlation_image_max(single(data.Y),8);

    % classify components
    [ROIvars.rval_space,ROIvars.rval_time,ROIvars.max_pr,ROIvars.sizeA,keep] = classify_components(data,A,C,b,f,YrA,options);

else
    % Run on patches (around 15 minutes)
    % keyboard
    [A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(Y,K,patches,tau,p,options);

    % compute correlation image on a small sample of the data (optional - for visualization purposes)
    Cn = correlation_image_max(single(Y),8);

    % classify components
    [ROIvars.rval_space,ROIvars.rval_time,ROIvars.max_pr,ROIvars.sizeA,keep] = classify_components(Y,A,C,b,f,YrA,options);

end

%% run GUI for modifying component selection (optional, close twice to save values)
Coor    = plot_contours(A,Cn,options,1); close;
keep    = filter_border_ROIs_jc(keep, FOV, Coor);
run_GUI = false;
if run_GUI
    Coor    = plot_contours_old(A,Cn,options,1); close;
    GUIout  = ROI_GUI_old(A,options,Cn,Coor,keep,ROIvars);
    options = GUIout{2};
    keep    = GUIout{3};
end

%% view contour plots of selected and rejected components (optional)
throw = ~keep;
figure;
ax1 = subplot(121); plot_contours_old(A(:,keep),Cn,options,1,[],Coor,1,find(keep)); title('Selected components','fontweight','bold','fontsize',14);
ax2 = subplot(122); plot_contours_old(A(:,throw),Cn,options,1,[],Coor,1,find(throw));title('Rejected components','fontweight','bold','fontsize',14);
linkaxes([ax1,ax2],'xy')

%% plot all rois including manually drawn ones
figure,
plot_contours_old(A,Cn,options,1,[],Coor,1,1:size(A,2)); title('Selected components','fontweight','bold','fontsize',14);

%% keep only the active components
A_keep = A(:,keep);
C_keep = C(keep,:);

%% extract fluorescence and DF/F on native temporal resolution
% C is deconvolved activity, C + YrA is non-deconvolved fluorescence
% F_df is the DF/F computed on the non-deconvolved fluorescence
files = dir(fullfile([sessionfoldername 'PreProcess' filesep area_channel filesep 'Avg' filesep ],['*' area_channel '*.tif']));   % list of Avg files
numFiles    = length(files);

P.p    = 0;                 % order of dynamics. Set P.p = 0 for no deconvolution at the moment
C_us   = cell(numFiles,1);  % cell array for thresholded fluorescence
f_us   = cell(numFiles,1);  % cell array for temporal background
P_us   = cell(numFiles,1);
S_us   = cell(numFiles,1);
YrA_us = cell(numFiles,1);  %
b_us   = cell(numFiles,1);  % cell array for spatial background

if ~exist('Ts','var')
    Ts     = data.Ts;        % added by Jianguang
end
parfor i = 1:numFiles
    
    fullname          = files(i).name;          % added by Jerry
    [~,file_name,ext] = fileparts(fullname);    % added by Jerry
    file_name         = file_name(5:end);
    folder_name       = [sessionfoldername 'PreProcess' filesep area_channel filesep ];
    name              = fullfile(folder_name,[file_name,ext]);
    
    [C_us{i},f_us{i},~,~,YrA_us{i}] = update_temporal_components_axg(name,A_keep,b,[],[],P,options);
    b_us{i} = max(mm_fun(f_us{i},name) - A_keep*(C_us{i}*f_us{i}'),0)/norm(f_us{i})^2;
end

prctfun = @(data) prctfilt(data,30,1000,300);       % first detrend fluorescence (remove 30th percentile on a rolling 1000 timestep window)
F_us    = cellfun(@plus,C_us,YrA_us,'un',0);        % cell array for projected fluorescence
Fd_us   = cellfun(prctfun,F_us,'un',0);             % detrended fluorescence

Ab_d = cell(numFiles,1);                            % now extract projected background fluorescence
for i = 1:numFiles
    Ab_d{i} = prctfilt((bsxfun(@times, A_keep, 1./sum(A_keep.^2))'*b_us{i})*f_us{i},30,1000,300,0);
end

F0 = cellfun(@plus, cellfun(@(x,y) x-y,F_us,Fd_us,'un',0), Ab_d,'un',0);   % add and get F0 fluorescence for each component
% DF/F value


%% export all necessary files for step 3
% Use F_dF and Coor
close all;

F_df    = cellfun(@(x,y) x./y, Fd_us, F0 ,'un',0);
ROIs    = Coor(keep);
% ROIs    = Coor;
[CC,im] = plot_contours2(Cn,1,ROIs,reference_image, sessionfoldername, ref_trial);

% /// generate cell labels
for i = 1:length(ROIs)
    cellid{i,1}     = [sessionfile(1:end-4) '-' area_channel(1:2) '-' sprintf('%04d', i) '-temp'];
    celltype{i,1}   = 'U';
    cellalgo{i,1}   = 'CNMF';       % added by Jianguang
end


%% save data

% /// check and save trial_info : added by Jianguang
files = dir(fullfile([sessionfoldername 'PreProcess' filesep area_channel filesep 'Avg' filesep],['*' area_channel '*.tif']));
if exist('trial_info','var')
    
    if length(trial_info) == length(files)
        Ca.trial_info = trial_info;
    end
    
else
    % /// read trial info from raw data
    for i = 1:length(files)
        filename   = files(i).name;
        filename   = filename(5:end-4);
        foldername = [sessionfoldername 'PreProcess' filesep area_channel filesep];
        s = load([foldername filename '.mat']);
        trial_info{i}.motion_metric = s.motion_metric;
        trial_info{i}.time_stamp    = s.time_stamp;
        trial_info{i}.fileloc       = s.fileloc;
        trial_info{i}.mat_file      = [filename '.mat'];
    end
    Ca.trial_info = trial_info;
    
end

% /// save other information
Ca.b          = b;
Ca.F_dF       = F_us;
Ca.ROIs       = ROIs;
Ca.trial_info = trial_info;
Ca.cellid     = cellid;
Ca.celltype   = celltype;
Ca.celalgo    = cellalgo;
Ca.ref_trial  = ref_trial;

if isempty(findstr(area_channel, 'A0')) == 0
    CaA0 = Ca;
    CaA0.sessionfoldername = strrep(Ca.sessionfoldername,'/net/claustrum2/mnt/data','X:');
    save([savefoldername sessionfile], 'CaA0','-append');
    disp(['successfully save CaA0 to ' sessionfile])
    
elseif isempty(findstr(area_channel, 'A1')) == 0
    CaA1 = Ca;
    CaA1.sessionfoldername = strrep(Ca.sessionfoldername,'/net/claustrum2/mnt/data','X:');
    save([savefoldername sessionfile], 'CaA1','-append');
    
elseif isempty(findstr(area_channel, 'A2')) == 0
    CaA2 = Ca;
    save([savefoldername sessionfile], 'CaA2','-append');
    disp(['successfully save CaA2 to ' sessionfile])
    
elseif isempty(findstr(area_channel, 'A3')) == 0
    CaA3 = Ca;
    save([savefoldername sessionfile], 'CaA3','-append');
end

delete(ds_filename);
poolobj = gcp('nocreate');
delete(poolobj);

if exist('slackID','var')
    if ~isempty(slackID)
        SendSlackNotification('https://hooks.slack.com/services/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', ...
            ['<@', slackID, '> 2P PIPELINE STEP 02 ' anm '-' sessionNo subsession ' ' area_channel ' is finished'], '#e_pipeline_log');
    else
        SendSlackNotification('https://hooks.slack.com/services/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', ...
            ['2P PIPELINE STEP 02 ' anm '-' sessionNo subsession ' ' area_channel ' is finished'], '#e_pipeline_log');
    end
else
    SendSlackNotification('https://hooks.slack.com/services/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', ...
        ['2P PIPELINE STEP 02 ' anm '-' sessionNo subsession ' ' area_channel ' is finished'], '#e_pipeline_log');
end
end

function [fileloc, mat_file, time_stamp, motion_metric] = trial_align_batch(fullname, sessionfoldername, area_channel, reference_image, redo_trial_alignment, multi_plane)

reg_image = imread([sessionfoldername 'PreProcess' filesep area_channel filesep 'Avg' filesep fullname]);
mat_file = [fullname(5:end-4) '.mat'];
load([sessionfoldername 'PreProcess' filesep area_channel filesep mat_file]);

if exist('aligned','var') == 0
   aligned = 0; 
end

if multi_plane
    num_planes=2;
    reference_image1=reference_image(1:end/2,:);
    reference_image2=reference_image(end/2+1:end,:);
    reg_image1=reg_image(1:end/2,:);
    reg_image2=reg_image(end/2+1:end,:);
    motion_corrected1=motion_corrected(1:end/2,:,:);
    motion_corrected2=motion_corrected(end/2+1:end,:,:);
else
    num_planes=1;
end

for i=1:num_planes
    if multi_plane
        if i==1
            offset1    = align_image_batch(reference_image1, reg_image1);
            if aligned == 0 || redo_trial_alignment == 1
                for j = 1:size(motion_corrected,3)
                    motion_corrected1(:,:,j) = circshift(motion_corrected1(:,:,j), [offset1(2) offset1(1)]);
                end
            end
            save([sessionfoldername 'PreProcess' filesep area_channel filesep mat_file], 'offset1', '-append');
        else
            offset2    = align_image_batch(reference_image2, reg_image2);
            if aligned == 0 || redo_trial_alignment == 1
                for j = 1:size(motion_corrected,3)
                    motion_corrected2(:,:,j) = circshift(motion_corrected2(:,:,j), [offset2(2) offset2(1)]);
                end
            end
            
            aligned = 1;
            motion_corrected=cat(1,motion_corrected1,motion_corrected2);
            save([sessionfoldername 'PreProcess' filesep area_channel filesep mat_file], 'motion_corrected', '-append');
            save([sessionfoldername 'PreProcess' filesep area_channel filesep mat_file], 'offset2', '-append');
            save([sessionfoldername 'PreProcess' filesep area_channel filesep mat_file], 'aligned', '-append');
        end
    else
        offset    = align_image_batch(reference_image, reg_image);
        
        if aligned == 0 || redo_trial_alignment == 1
            for j = 1:size(motion_corrected,3)
                motion_corrected(:,:,j) = circshift(motion_corrected(:,:,j), [offset(2) offset(1)]);
            end
            
            aligned = 1;
            save([sessionfoldername 'PreProcess' filesep area_channel filesep mat_file], 'motion_corrected', '-append');
            save([sessionfoldername 'PreProcess' filesep area_channel filesep mat_file], 'offset', '-append');
            save([sessionfoldername 'PreProcess' filesep area_channel filesep mat_file], 'aligned', '-append');
        end
    end
end


tiff_info.directory     = [sessionfoldername 'PreProcess' filesep area_channel filesep];
tiff_info.filename      = [fullname(5:end-4) '.tif'];
tiff_info.bitspersample = 16;
tiff_info.floatingpoint = false;
tiff_info.nframes       = size(motion_corrected, 3);

UCLA_WriteTiff(uint16(motion_corrected), tiff_info);
end
