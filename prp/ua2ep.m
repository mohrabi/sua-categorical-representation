close all
clear

monkey = "Jenab";
selectivity = "Fast";

monkey_dir = fullfile('G:\Data', selectivity, monkey);
session_dates = ls(fullfile(monkey_dir, '20*'));

f = fopen('G:\Codes\Preprocessing\Log\ua2ep_log.txt', 'at+');
fprintf(f, strcat(datestr(datetime('now')), '\n'));
fprintf(f, "\tGenerating data epochs from sorting results...\n");

global pre pos
if selectivity == "Fast"
    pre = 200;
    pos = 699;
elseif selectivity == "Slow"
    pre = 200;
    pos = 999;
end

%%
tic
for i_session = 113:size(session_dates, 1)
    clc, disp("Creating data epochs for monkey " + monkey + ...
        ", session " + num2str(i_session) + "/" + num2str(size(session_dates, 1)))
    fprintf(f, "\tSession <" + session_dates(i_session, :) + ...
        "> of monkey <" + monkey + ">\n");
    
    try
        clean_dir = fullfile(monkey_dir, session_dates(i_session, :), 'Sort');
        trials_dir = fullfile(monkey_dir, session_dates(i_session, :), 'Trial');
        mat_dir = fullfile(monkey_dir, session_dates(i_session, :), 'Mat');
        
        if isfolder(trials_dir), rmdir(trials_dir, 's'), end
        mkdir(trials_dir)
        
        load(fullfile(mat_dir, 'evt.mat'))
        [trial_time, trial_info] = parse_events(tim, evt);
        clear tim evt
        
        clear cm
        cm = create_condition_matrix(trial_info);
        save(fullfile(trials_dir, "cm.mat"), 'cm');
                
        if ~isfile(fullfile(clean_dir, '__IT_NO_LOAD__'))
            [it_spike_times, it_optimal_set] = load_sorted_spikes(fullfile(clean_dir, 'it.mat'));
            it_epochs = create_data_epochs(it_spike_times, it_optimal_set, trial_time);
            save_epochs(it_epochs, trials_dir, 'it')
        else
            fprintf(f, "\t\tSkipping IT data.\n");
        end
        clear it_spike_times it_optimal_set it_epochs
        
        if ~isfile(fullfile(clean_dir, '__PFC_NO_LOAD__'))
            [pfc_spike_times, pfc_optimal_set] = load_sorted_spikes(fullfile(clean_dir, 'pfc.mat'));
            pfc_epochs = create_data_epochs(pfc_spike_times, pfc_optimal_set, trial_time);
            save_epochs(pfc_epochs, trials_dir, 'pfc')
        end
        clear pfc_spike_times pfc_optimal_set pfc_epochs
        
        % if ~isfile(fullfile(clean_dir, '__PFC_NO_LOAD__'))
        %     if isfile(fullfile(clean_dir, 'pfc.mat'))
        %         [pfc_spike_times, pfc_optimal_set] = load_sorted_spikes(fullfile(clean_dir, 'pfc.mat'));
        %     else
        %         [pfc_spike_times, pfc_optimal_set] = load_sorted_spikes(fullfile(mat_dir, 'spk_pfc_sor.mat'));
        %     end
        %     pfc_epochs = create_data_epochs(pfc_spike_times, pfc_optimal_set, trial_time);
        %     save_epochs(pfc_epochs, trials_dir, 'pfc')
        % end
        
        
        if length(cm) < 50
            fprintf(f, "\t\tSession <" + session_dates(i_session, :) + ...
                "> of monkey <" + monkey + "> has only " + length(cm) + " trials.\n");
        end
        
    catch ME
        fprintf(f, "\t\t" + ME.message + ".\n");
    end
end
fprintf(f, "\tFinished successfully after " + num2str(round(toc)) + ...
    " seconds.\n\n");
fclose(f);

%%
% i_session = 60;
% clean_dir = fullfile(monkey_dir, session_dates(i_session, :), 'Sort');
% trials_dir = fullfile(monkey_dir, session_dates(i_session, :), 'Trial');
% mat_dir = fullfile(monkey_dir, session_dates(i_session, :), 'Mat');

%%
% if isfolder(trials_dir), rmdir(trials_dir, 's'), end
% mkdir(trials_dir)
% 
% load(fullfile(mat_dir, 'evt.mat'))
% [trial_time, trial_info] = parse_events(tim, evt);
% clear tim evt
%%
% clear cm
% cm = create_condition_matrix(trial_info);
% save(fullfile(trials_dir, "cm.mat"), 'cm');
%%
% if ~isfile(fullfile(clean_dir, '__IT_NO_LOAD__'))
%     [it_spike_times, it_optimal_set] = load_sorted_spikes(fullfile(clean_dir, 'it.mat'));
%     it_epochs = create_data_epochs(it_spike_times, it_optimal_set, trial_time);
%     save_epochs(it_epochs, trials_dir, 'it')
% else
%     fprintf(f, "\t\tSkipping IT data.\n");
% end
%%
% if ~isfile(fullfile(clean_dir, '__PFC_NO_LOAD__'))
%     [pfc_spike_times, pfc_optimal_set] = load_sorted_spikes(fullfile(clean_dir, 'pfc.mat'));
%     pfc_epochs = create_data_epochs(pfc_spike_times, pfc_optimal_set, trial_time);
%     save_epochs(pfc_epochs, trials_dir, 'pfc')
% end
%%
% if length(cm) < 50
%     fprintf(f, "\t\tSession <" + session_dates(i_session, :) + ...
%         "> of monkey <" + monkey + "> has only " + length(cm) + " trials.\n");
% end
    
%%


function [spike_times, optimal_set] = load_sorted_spikes(file_name)
spike_times  = load(file_name);
optimal_set = spike_times.optimal_set.cluster_index;
spike_times  = round(spike_times.SpikeTime / 40);
end

function [event_time, event_tag] = parse_events(tim, evt)
tim = tim + tim(1);
tim = tim(3:end);
evt = evt(3:end);

ix_h = (evt==0 | evt==255);
evt(ix_h) = [];
tim(ix_h) = [];

ix = find(ismember(evt, 1:170));

event_tag = []; event_time = [];
for t = 1:(length(ix)-2)
    if evt(ix(t))+55 == evt(ix(t)+1) && evt(ix(t))+75 == evt(ix(t)+2)
        event_tag = [event_tag; evt(ix(t))]; %#ok<AGROW>
        event_time = [event_time; tim(ix(t))]; %#ok<AGROW>
    end
end
event_time = round(event_time / 40);
end

function epochs = create_data_epochs(spike_times, optimal_set, onsets)
global pre pos

spike_times(optimal_set > 250) = [];
optimal_set(optimal_set > 250) = [];
optimal_set_ids = unique(optimal_set);
nu = length(optimal_set_ids);
last_spike_sample = max(spike_times);
epochs = false(size(onsets, 1), nu, pre + pos + 1);

for iu = 1:nu
    spike_train = false(last_spike_sample + 1, 1);
    spike_train(spike_times(optimal_set == iu)) = true;
    
    for i_trial = 1:length(onsets)
        [start_time, finish_time] = get_trial_time_markers(onsets(i_trial));
        
        if start_time < 0 || finish_time > onsets(end), continue, end
        if(finish_time > last_spike_sample), break, end
        epochs(i_trial, iu, :) = spike_train(start_time:finish_time);
    end
    if(i_trial ~= length(onsets))
        epochs(i_trial:end, :, :) = [];
    end
end
end

function [tar, fin] = get_trial_time_markers(onset)
global pre pos
tar  = onset - pre;
fin = onset + pos;
end

function cm = create_condition_matrix(trial_info)
cm = trial_info;

% cm = nan(size(trialInf, 1), 11); % ind-art-nat-npf-mof-huf-hub-anb-bod-fac-ina-ani
% ART = 2; NAT = 3; NPF = 4; MOF = 5; HUF = 6; HUB = 7; ANB = 8; BOD = 9;
% FAC = 10; INA = 11; ANI = 12; LSF = 13; ISF = 14; HSF = 15;

% for iTrial = 1:length(trialTim)
%     cm(iTrial, 1)   = trialInf(iTrial);
%     c               = trialInf(iTrial);
%     cm(iTrial, ART) = any(c == [55:71,74,96:101,123:128]);
%     cm(iTrial, NAT) = any(c == [38:54,72,73,90:95,116:122,144:149]);
%     cm(iTrial, NPF) = any(c == [10:12,14:16,81,82,108,109,135,136]);
%     cm(iTrial, MOF) = any(c == [13,17,18,83,110,137,166:170]);
%     cm(iTrial, HUF) = any(c == [1:9,75:80,102:107,129:134,156:165]);
%     cm(iTrial, HUB) = any(c == [19:28,84:86,111:113,138:140]);
%     cm(iTrial, ANB) = any(c == [29:37,87:89,114:116,141:143]);
%     cm(iTrial, BOD) = any(cm(iTrial, [ANB,HUB]));
%     cm(iTrial, FAC) = any(cm(iTrial, [NPF,MOF,HUF]));
%     cm(iTrial, INA) = any(cm(iTrial, [ART,NAT]));
%     cm(iTrial, ANI) = any(cm(iTrial, [FAC,BOD]));
%     cm(iTrial, HSF) = any(c == 75:101);
%     cm(iTrial, LSF) = any(c == 129:155);
%     cm(iTrial, ISF) = ~(cm(iTrial, HSF) || cm(iTrial, LSF));
% end
% clear ART NAT NPF MOF HUF HUB ANB BOD FAC INA ANI LSF HSF ISF
% cm = array2table(cm, 'VariableNames', ...
%         ["ind", "art", "nat", "npf", "mof", "huf" ,"hub", "anb", "bod", "fac", "ina", "ani", "lsf", "isf", "hsf"]);
end

function save_epochs(epochs, save_dir, channel_name)
n_spike = sum(epochs, [1, 3]);
[~, mu_index] = max(n_spike);
ua = squeeze(epochs(:, mu_index, :)); %#ok<*NASGU>
file_name = strcat(strjoin({'mu', channel_name}, '_'), '.mat');
save(fullfile(save_dir, file_name), 'ua')
epochs(:, mu_index, :) = [];

if isempty(epochs), return, end
if ndims(epochs) == 3
    for iu = 1:size(epochs, 2)
        ua = squeeze(epochs(:, iu, :));
        file_name = strcat(strjoin({'su', channel_name, num2str(iu)}, '_'), '.mat');
        save(fullfile(save_dir, file_name), 'ua')
    end
else
    error('save_epochs: Incorrect matrix shape')
end

end