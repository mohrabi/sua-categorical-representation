close all
clear
clc

selectivity = "Fast";

monkey_dir = fullfile('J:\Data', selectivity, 'Jenab');
export_dir = fullfile('G:\Data', selectivity, 'Jenab');

sessions = readmatrix('G:\Data\Fast\Jenab\recordingInfo.csv', ...
    'OutputType', 'string', ...
    'NumHeaderLines', 1, ...
    'ExpectedNumVariables', 2, ...
    'Delimiter', ',');
sessions = sessions(sessions(:, 2) == "TRUE", 1);
sessions = replace(sessions, "'", "");
sessions = string(datestr(...
    datetime(sessions, 'InputFormat', "dd-MMM-yyyy HH:mm:ss"), ...
    "yyyy-mm-dd_HH-MM"));

nsession = size(sessions, 1);

global fs pre pos
fs = 1000;
if selectivity == "Fast"
    pre = 1200;
    pos = 1699;
elseif selectivity == "Slow"
    pre = 1200;
    pos = 1999;
end

for isession = 1:nsession
    clc, disp([isession, nsession])
    mat_dir = fullfile(monkey_dir, sessions(isession, :), 'mat');
    lfp_dir = fullfile(export_dir, sessions(isession, :), 'Trial');
    
    % Load and Save Condition Matrix
    load(fullfile(mat_dir, 'evt.mat'))
    [trial_time, ~] = parse_events(tim, evt);
    clear tim evt
    
    lfp = load_local_field_potentials(fullfile(mat_dir, 'lfp_it.mat'));
    data = create_data_epochs(lfp, trial_time);
    save(fullfile(lfp_dir, 'l_it.mat'), 'data')
    
    lfp = load_local_field_potentials(fullfile(mat_dir, 'lfp_pfc.mat'));
    data = create_data_epochs(lfp, trial_time);
    save(fullfile(lfp_dir, 'l_pfc.mat'), 'data')
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

function signal = load_local_field_potentials(file_name)
signal  = load(file_name);
signal = struct2array(signal);
signal = decimate(signal, 5);
signal = decimate(signal, 8);

load bpdf.mat
signal = filtfilt(bpdf, signal);

end

function [tar, fin] = get_trial_time_markers(onset)
global pre pos fs
tar  = onset - pre / 1000 * fs;
fin = onset + pos / 1000 * fs;
end

function epochs = create_data_epochs(lfp, onsets)
global pre pos fs
epochs = NaN(size(onsets, 1), (pre + pos) * fs / 1000 + 1);
for i_trial = 1:length(onsets)
    [start_time, finish_time] = get_trial_time_markers(onsets(i_trial));
    
    if start_time < 0 || finish_time > onsets(end), break, end
    if finish_time > size(lfp, 1), break, end
    epochs(i_trial, :) = lfp(start_time:finish_time);
end
if(i_trial ~= length(onsets))
    epochs(i_trial:end, :) = [];
end
end
