close all
clear
clc

addpath('G:\Codes\Common')
selectivity = "Slow";
monkey = "Zebel";

root = fullfile('G:\Data', selectivity, monkey);

sessions = ls(fullfile(root, '20*'));
nsession = size(sessions, 1);
%%
info = ["file path", "unit id", ...
    "overal p-value", "Face p-value", "Body p-value", ...
    "Artificial p-value", "Natural p-value"];
t = -200:699;

for channel = ["it", "pfc"]
    for isession = 1:nsession
        trial_dir = fullfile(root, sessions(isession, :), 'Trial');
        load(fullfile(trial_dir, "cm.mat"))
        
        try
            if isfile(fullfile(trial_dir, "mu_" + channel + ".mat"))
                clear ua, load(fullfile(trial_dir, "mu_" + channel + ".mat"))
                info = [info; fullfile(trial_dir, "mu_" + channel + ".mat"), ...
                    strcat(upper(channel), "/", "mu", "/", num2str(size(info, 1))), ...
                    per_category_significance(ua, cm(1:size(ua, 1)))];
            else
                error("epoched data file path is not valid!")
            end
        catch ME
            error(ME.message)
        end
        
        singleunits = ls(fullfile(trial_dir, "su_" + channel + "*" + ".mat"));
        nsingleunit = size(singleunits, 1);
        for isingleunit = 1:nsingleunit
            clear ua, load(fullfile(trial_dir, singleunits(isingleunit, :)))
            info = [info; fullfile(trial_dir, singleunits(isingleunit, :)), ...
                    strcat(upper(channel), "/", "su", "/", num2str(size(info, 1))), ...
                    per_category_significance(ua, cm(1:size(ua, 1)))];
        end
        clear singleunits nsingleunit isingleunit
    end
end
clear channel isession isingleunit
clear ua

%% Dont need the p-values
info = info(:, 1:2);
%% Save but not overwrite
if isfile(fullfile(root, "guide.csv"))
    error("file exists")
else
    writematrix(info, fullfile(root, "guide.csv"))
end
%%
% %%
% recording_info = readmatrix(fullfile('G:\Data', "Fast", "Jenab", "Jenab.xlsx"), ...
%     'FileType', 'spreadsheet', ...
%     'NumHeaderLines', 1, ...
%     'OutputType', 'string');
% foldername = recording_info(:, 1);
% dv_it = str2double(recording_info(:, 3));
% dv_pfc = str2double(recording_info(:, 2));
% 
% info_it = info(info(:, 1) == "IT", :);
% info_pfc = info(info(:, 1) == "PFC", :);
% 
% %%
% invalid_index = false(size(info_it, 1), 1);
% for i = 2:size(foldername, 1)
%     if dv_it(i) == dv_it(i-1)
%         invalid_index(datetime(info_it(:, 2), 'InputFormat', 'yyyy-MM-dd_HH-mm') == ...
%             datetime(extractBetween(foldername(i), 2, 21), 'InputFormat', "dd-MMM-yyyy HH:mm:ss")) = true;
%     end
% end
% clear i
% sessioninfo_it = info_it(~invalid_index, :);
% 
% invalid_index = false(size(info_pfc, 1), 1);
% for i = 2:size(foldername, 1)
%     if dv_pfc(i) == dv_pfc(i-1)
%         invalid_index(datetime(info_pfc(:, 2), 'InputFormat', 'yyyy-MM-dd_HH-mm') == ...
%             datetime(extractBetween(foldername(i), 2, 21), 'InputFormat', "dd-MMM-yyyy HH:mm:ss")) = true;
%     end
% end
% clear i
% sessioninfo_pfc = info_pfc(~invalid_index, :);
% 
% sessioninfo = [sessioninfo_it; sessioninfo_pfc];
% clear sessioninfo_it sessioninfo_pfc info_it info_pfc
% %%
% sessioninfo(all(ismissing(sessioninfo(:, end-4:end)), 2), :) = [];
% %%
% sessioninfo = ["region", "folder name", "folder index", "unit type", ...
%     "overal p-value", "Face p-value", "Body p-value", ...
%     "Artificial p-value", "Natural p-value"; sessioninfo];
% %%
% save(fullfile('G:\Data\Fast\Jenab\sessioninfo.mat'), 'sessioninfo')
% writematrix(sessioninfo, fullfile('G:\Data\Fast\Jenab\sessioninfo.csv'))
% %% Remove non-selective IT neurons
% channel = "it";
% 
% % opval = str2double(sessioninfo(:, 5));
% fpval = str2double(sessioninfo(:, 6));
% bpval = str2double(sessioninfo(:, 7));
% apval = str2double(sessioninfo(:, 8));
% npval = str2double(sessioninfo(:, 9));
% 
% conf = 0.05;
% 
% indit = sessioninfo(:, 1) == upper(channel) & (fpval < conf | bpval < conf | ...
%     npval < conf | apval < conf);
% indit = ~indit & sessioninfo(:, 1) == upper(channel);
% nit = sum(sessioninfo(:, 1) == upper(channel));
% sessioninfo(indit, :) = [];
% % sessioninfo(indit, 1) = "Yep";
% 
% clear opval fpval bpval apval npval
% 
% disp(num2str(sum(indit)) + " units from total " + ...
%     num2str(nit) + " units of " + channel + ...
%     " were statistially indifferent from baseline.")
% %%
% channel = "pfc";
% 
% % opval = str2double(sessioninfo(:, 5));
% fpval = str2double(sessioninfo(:, 6));
% bpval = str2double(sessioninfo(:, 7));
% apval = str2double(sessioninfo(:, 8));
% npval = str2double(sessioninfo(:, 9));
% 
% indpfc = sessioninfo(:, 1) == upper(channel) & (fpval < 0.01 | bpval < 0.01 | ...
%     npval < 0.01 | apval < 0.01);
% indpfc = ~indpfc & sessioninfo(:, 1) == upper(channel);
% npfc = sum(sessioninfo(:, 1) == upper(channel));
% sessioninfo(indpfc, :) = [];
% 
% clear opval fpval bpval apval npval
% 
% disp(num2str(sum(indpfc)) + " units from total " + ...
%     num2str(npfc) + " units of " + channel + ...
%     " were statistially indifferent from baseline.")
% %%
% save(fullfile('G:\Data\Fast\Jenab\selective_sessioninfo.mat'), 'sessioninfo')
% writematrix(sessioninfo, fullfile('G:\Data\Fast\Jenab\selective_sessioninfo.csv'))
% %%
function sig = per_category_significance(ua, cm)
t = -200:699;
ind1 = (t <= 30) & (t > -70);
ind2 = (t <= 170) & (t > 70);

sig = NaN(1, 5); % overall, face, body, artificial, natural
[sig(1), ~] = signrank(mean(ua(:, ind1), 2), mean(ua(:, ind2), 2));

y = grablabels('face-body');
z = ismember(cm, find(y == 0));
[sig(2), ~] = signrank(mean(ua(z, ind1), 2), mean(ua(z, ind2), 2));
z = ismember(cm, find(y == 1));
[sig(3), ~] = signrank(mean(ua(z, ind1), 2), mean(ua(z, ind2), 2));

y = grablabels('artificial-natural');
z = ismember(cm, find(y == 0));
[sig(4), ~] = signrank(mean(ua(z, ind1), 2), mean(ua(z, ind2), 2));
z = ismember(cm, find(y == 1));
[sig(5), ~] = signrank(mean(ua(z, ind1), 2), mean(ua(z, ind2), 2));
end