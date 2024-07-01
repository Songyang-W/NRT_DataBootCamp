function pipeline_Step2_SCC(anm,sessionList)
% pipeline_Step2.m
% Songyang
% Runs CHEN_2P_Pipeline_Step_02 for the specified animal and sessions
% Copy from David's pr_pipeline_Step2.m

TODO
%% Fill in the path directory (ending at /Animals/), animal name, and session list
pathdirectory = '/net/claustrum2/mnt/data/Projects/Connectomics/Animals/';
%anm = 'jc105';
%sessionList = 1;
SCCdir = pwd;
%% Do not make edits beneath this line

area_channels  = {'A0_Ch0', 'A1_Ch0'};
sampling_rate = 32.5868;                         % Hz, full frame acqusition rate of the resonance scanner
FOV           = [400 460];                       % pixel, field of view
do_parallel   = 1;                               % parallel computing option
do_nonrigid   = 1;                               % option to do non-rigid registration (faster without it, better results with it)
multi_plane   = 0;
downsample_f = 10;
redo_trial_alignment = 1;

pd = [pathdirectory anm '/'];


if isstr(sessionList)
    sessionList = str2num(sessionList);
end

pause(60 + randi([-30 30],1))

load([pd 'refs.mat'])
for sN = sessionList
    for aC = 1:2
        close all
        sessionNo = num2str(sN);
        sessionName   = [anm , '-', sessionNo];
        area_channel = area_channels{aC};
        aNum = area_channel(1:2);
        
        % Make sure reflist is up to date before checking/making changes
        clear('reflist','autoRef')
        load([pd 'refs.mat'])
        
        % Run Step 2 if it is labeled 'Ready' (manually Excluded Trials)
        if strcmp(reflist(sN).(aNum),'Ready')
            % Use assigned ref trial if it exists; else find automatically
            if ~isempty(reflist(sN).([aNum 'ref']))
                ref_trial = reflist(sN).([aNum 'ref']);
            else
                ref_trial = '';
            end

            reflist(sN).(aNum) = 'In Progress';
            save([pd 'refs.mat'],'reflist')
            
            % Run Step 2
            cd(SCCdir)
            CHEN_2P_Pipeline_Step_02(pathdirectory,anm, sessionNo, area_channel, ref_trial, ...
                sampling_rate, downsample_f, redo_trial_alignment, multi_plane,'')
            
            % Make sure reflist is up to date before making changes
            clear('reflist')
            
            % Find reference trial used if ref_trial is not specified
            if isempty(ref_trial)
                cd([pathdirectory anm '/2P/' anm '-' num2str(sN) '/PreProcess'])
                fList = ls;
                [numFs, ~] = size(fList);
                % Scan through files for Area # and .tif
                for fN = 1:numFs
                    if contains(fList(fN,:),'tif') && contains(fList(fN,:),area_channel) && ~contains(fList(fN,:),'activity')
                        CNMF_file = fList(fN,:);
                        autoRef = ['Avg_' CNMF_file(6:end)]; % Reference is listed as CNMF_A#[...].tif so remove the first 5 characters
                    end
                end
            end
            load([pd 'refs.mat'])
            
            % Update reflist
            reflist(sN).(aNum) = 'X';
            if exist('autoRef','var')
                reflist(sN).([aNum 'ref']) = autoRef;
            end
            % Save reflist in original location
            save([pd 'refs.mat'],'reflist')
        end
    end
end
end
