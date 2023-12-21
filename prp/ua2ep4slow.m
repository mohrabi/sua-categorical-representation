close all
clear

monkey = "Zebel";
selectivity = "Slow";

monkey_dir = fullfile('J:\Data', selectivity, monkey);
export_dir = fullfile('G:\Data', selectivity, monkey);
sessions = ls(fullfile(monkey_dir, '20*'));
nsession = size(sessions, 1);

global pre pos
if selectivity == "Fast"
    pre = 200;
    pos = 699;
elseif selectivity == "Slow"
    pre = 200;
    pos = 999;
end

%%
for isession = 1:nsession
    disp([isession, nsession])
    
    % Define Directories
    cln_dir = fullfile(monkey_dir, sessions(isession, :), 'sorted');
    mat_dir = fullfile(monkey_dir, sessions(isession, :), 'mat');
    
    % Create Outputs' Directory
    out_dir = fullfile(export_dir, sessions(isession, :), 'Trial');
    if ~isfolder(out_dir)
        mkdir(out_dir)
    end
    
    % Load and Save Condition Matrix
    load(fullfile(mat_dir, 'evt.mat'))
    [trial_time, trial_info] = parse_events(tim, evt);
    clear tim evt cm
    cm = create_condition_matrix(trial_info);
    save(fullfile(out_dir, "cm.mat"), 'cm');
    
    if length(cm) < 100
        disp(sessions(isession, :))
    end
    
    % Load IT Trials
    [it_spike_times, it_optimal_set] = load_sorted_spikes(fullfile(cln_dir, 'it.mat'));
    it_epochs = create_data_epochs(it_spike_times, it_optimal_set, trial_time);
    save_epochs(it_epochs, out_dir, 'it')
    
    % Load PFC Trials
    [it_spike_times, it_optimal_set] = load_sorted_spikes(fullfile(cln_dir, 'pfc.mat'));
    it_epochs = create_data_epochs(it_spike_times, it_optimal_set, trial_time);
    save_epochs(it_epochs, out_dir, 'pfc')
end

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

function epochs = create_data_epochs_from_mu(spike_times, optimal_set, onsets)
global pre pos

% load only the whole multi unit
spike_times(optimal_set > 200) = [];
optimal_set(optimal_set > 200) = [];
optimal_set(:) = 1;

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