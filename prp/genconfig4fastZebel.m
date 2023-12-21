%% Load Recording Info

close all
clear
clc

addpath('G:\Codes\Common')
addpath('C:toolbox\Separation Toolbox')

root_good = fullfile('J:\Data', "Fast", "Zebel");
root_corr = fullfile('J:\Data', "Fast", "Zebel", 'Corrupted');

sessions = [string(ls(fullfile(root_good, '20*')));
    string(ls(fullfile(root_corr, '20*')))];
sessions = sort(sessions);
nsession = size(sessions, 1);

sessionInfo = sessions;
sessionInfo(:, 2) = "N/A";
sessionInfo(ismember(sessions, string(ls(fullfile(root_good, '20*')))), 2) = "good";
sessionInfo(ismember(sessions, string(ls(fullfile(root_corr, '20*')))), 2) = "corrupted";

recordingInfo = readmatrix(fullfile('G:\Data', "Fast", "Zebel", "recordingInfo.xlsx"), ...
    'FileType', 'spreadsheet', ...
    'NumHeaderLines', 0, ...
    'OutputType', 'string');
if size(recordingInfo, 2) < 8
    recordingInfo(:, 4:8) = missing;
elseif size(recordingInfo, 2) < 8
    error("")
end

dv = str2double(recordingInfo(:, 2:3));
hole = recordingInfo(:, [4, 8]);
sessionInfo(:, 3:6) = [dv, hole];

clear dv hole recordingInfo sessions
writematrix(sessionInfo, 'G:\Data\Fast\Zebel\sessionInfo.csv');

%% Generate Neuron-Specific Info
clear
clc
close all

sessionInfo = readmatrix('G:\Data\Fast\Zebel\sessionInfo.csv', ...
    'FileType', 'spreadsheet', ...
    'NumHeaderLines', 0, ...
    'OutputType', 'string');
info = sessionInfo(sessionInfo(:, 2) ~= "corrupted", :);
if size(info, 2) < 6
    info(:, size(info, 2):6) = missing;
elseif size(info, 2) < 8
    error("")
end

itcNeuralInfo = [];
pfcNeuralInfo = [];
for iRow = 1:size(info, 1)
    root = fullfile('G:\Data\Fast\Zebel', info(iRow, 1), 'Trial');
    cm = load(fullfile(root, 'cm.mat'), 'cm').cm;
    itcNeurons = string(ls(fullfile(root, '*u*it*mat')));
    if itcNeurons ~= ""
        for iNeuron = size(itcNeurons, 1)
            ua = load(fullfile(root, itcNeurons(iNeuron)), 'ua').ua;
            sig = per_category_significance(ua, cm);
            if any(sig < 0.05)
                itcNeuralInfo = [
                    itcNeuralInfo; ...
                    fullfile(root, itcNeurons(iNeuron)), ...
                    extractBefore(itcNeurons(iNeuron), 3), ...
                    info(iRow, 2:end)
                    ]; %#ok<AGROW>
            end
        end
    end
    
    pfcNeurons = string(ls(fullfile(root, '*u*pfc*mat')));
    if pfcNeurons ~= ""
        for iNeuron = size(pfcNeurons, 1)
            ua = load(fullfile(root, pfcNeurons(iNeuron)), 'ua').ua;
            sig = per_category_significance(ua, cm);
            if any(sig < 0.05)
                pfcNeuralInfo = [
                    pfcNeuralInfo; ...
                    fullfile(root, pfcNeurons(iNeuron)), ...
                    extractBefore(pfcNeurons(iNeuron), 3), ...
                    info(iRow, 2:end)
                    ]; %#ok<AGROW>
            end
        end
    end
end

itcNeuralInfo = ["path", "type", "status", "dv1", "dv2", "h1", "h2"; itcNeuralInfo];
pfcNeuralInfo = ["path", "type", "status", "dv1", "dv2", "h1", "h2"; pfcNeuralInfo];

writematrix(itcNeuralInfo, 'G:\Data\Fast\Zebel\itcNeuralInfo.csv');
save('G:\Data\Fast\Zebel\itcNeuralInfo.mat', 'itcNeuralInfo');
writematrix(pfcNeuralInfo, 'G:\Data\Fast\Zebel\pfcNeuralInfo.csv');
save('G:\Data\Fast\Zebel\pfcNeuralInfo.mat', 'pfcNeuralInfo');

%%
figure_path = 'G:\Codes\Figures\psth\jenab-fast';
for iRow = 1:size(itcNeuralInfo, 1)
    ua = load(itcNeuralInfo(iRow, 1), 'ua').ua;
    cm = load(fullfile(fileparts(itcNeuralInfo(iRow, 1)), 'cm.mat'), 'cm').cm;
    
    myplot(-200:699, ua, cm)
    title("Fast Selectivity PSTH")
    subtitle("IT, ID: " + num2str(iRow))
    touch(gca)
    clear ua
    
    saveas(gcf, fullfile(figure_path, "IT-" + num2str(iRow) + ".png"))
    close gcf
end

%%
for iRow = 1:size(pfcNeuralInfo, 1)
    ua = load(pfcNeuralInfo(iRow, 1), 'ua').ua;
    cm = load(fullfile(fileparts(pfcNeuralInfo(iRow, 1)), 'cm.mat'), 'cm').cm;
    
    myplot(-200:699, ua, cm)
    title("Fast Selectivity PSTH")
    subtitle("PFC, ID: " + num2str(iRow))
    touch(gca)
    clear ua
    
    saveas(gcf, fullfile(figure_path, "pfc-" + num2str(iRow) + ".png"))
    close gcf
end

%%
% %%
% 
% info = ["region", "folder name", "folder index", "unit type", ...
%     "overal p-value", "Face p-value", "Body p-value", ...
%     "Artificial p-value", "Natural p-value"];
% t = -200:699;
% 
% for channel = ["it", "pfc"]
%     for isession = 1:nsession
%         trial_dir = fullfile(root, sessions(isession, :), 'Trial');
%         load(fullfile(trial_dir, "cm.mat"))
%         
%         try
%             if isfile(fullfile(trial_dir, "mu_" + channel + ".mat"))
%                 clear ua, load(fullfile(trial_dir, "mu_" + channel + ".mat"))
%                 info = [info; upper(channel), sessions(isession, :), isession, ...
%                     0, per_category_significance(ua, cm(1:size(ua, 1)))];
%             else
%                 info = [info; upper(channel), sessions(isession, :), isession, ...
%                     0, NaN, NaN, NaN, NaN, NaN];
%             end
%         catch ME
%             info = [info; upper(channel), sessions(isession, :), isession, ...
%                     0, NaN, NaN, NaN, NaN, NaN];
%         end
%         
%         singleunits = ls(fullfile(trial_dir, "su_" + channel + "*" + ".mat"));
%         nsingleunit = size(singleunits, 1);
%         for isingleunit = 1:nsingleunit
%             clear ua, load(fullfile(trial_dir, singleunits(isingleunit, :)))
%             info = [info; upper(channel), sessions(isession, :), isession, ...
%                 isingleunit, per_category_significance(ua, cm)];
%         end
%         clear singleunits nsingleunit isingleunit
%     end
% end
% clear channel isession isingleunit
% clear ua
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
%%
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

function myplot(t, ua, cm)
nanind = all(isnan(ua), 2);
figure('Units', 'centimeters', 'Position', [0, 0, 36, 15])
ax = subplot(1,1,1); hold on

mu = mean(ua(:, (t < 30) & (t > -70)), 'all');
sigma = std(ua(:, (t < 30) & (t > -70)), [], 'all');
ua = (ua - mu) ./ sigma;

y = grablabels('artificial-natural');
z0 = ismember(cm, find(y == 0));
z1 = ismember(cm, find(y == 1));
X = [ua(z0, :); ua(z1, :)];
groups = [zeros(sum(z0), 1); ones(sum(z1), 1)];
psth(ax, t, X, groups)


y = grablabels('face-body');
z0 = ismember(cm, find(y == 0));
z1 = ismember(cm, find(y == 1));
X = [ua(z0, :); ua(z1, :)];
groups = [zeros(sum(z0), 1); ones(sum(z1), 1)];
psth(ax, t, X, groups)

xlabel('Time(ms)')
ylabel('Fire Rate (Baseline Normalized)')
xline(0, 'k--', 'HandleVisibility', 'off', 'LineWidth', 2)
xlim([min(t), max(t)])

legend("Artificial", "Natural", "Face", "Body")
end

function psth(ax, t, X, groups)
group_numbers = unique(groups);
ngroup = length(group_numbers);
t = t(:);

X = smoothdata(X, 2, 'gaussian', 40);

for igroup = 1:ngroup
    alpha = .05; % Confidence Level
    g = X(groups == group_numbers(igroup), :);
    mu_g = mean(g, 1);
    s_g = std(g, [], 1);
    serr_g = tinv(alpha/2, size(g, 1)) * s_g / sqrt(size(g, 1));
    serr_g = serr_g(:); mu_g = mu_g(:);

    color = ax.ColorOrder(ax.ColorOrderIndex, :);
    plot(ax, t, mu_g , 'LineWidth', 3, 'Color', color);
    fill([t; flip(t)], [mu_g + serr_g; flip(mu_g-serr_g)], ...
        color, 'LineWidth', 3, 'FaceAlpha', .2, ...
        'EdgeColor', 'none', 'HandleVisibility', 'off')
    ax.ColorOrderIndex = mod(ax.ColorOrderIndex-2, size(ax.ColorOrder, 1))+1;
    
end
end