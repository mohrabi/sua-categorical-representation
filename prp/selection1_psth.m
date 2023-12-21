close all
clear
clc

monkey = "Zebel";
selectivity = "Slow";

monkey_dir = fullfile('G:\Data', selectivity, monkey);
figure_dir = 'G:\Codes\Figures\psth\zebel-slow';
sessions = ls(fullfile(monkey_dir, '20*'));
nsession = size(sessions, 1);

addpath('C:\toolbox\Separation Toolbox')
addpath('G:\Codes\Common')

if selectivity == "Fast"
    t = -200:699;
elseif selectivity == "Slow"
    t = -200:999;
end

%%
for channel = ["it", "pfc"]
    for isession = 1:nsession
        folder_name = datetime(sessions(isession, :), 'InputFormat', 'yyyy-MM-dd_HH-mm');
        trial_dir = fullfile(monkey_dir, sessions(isession, :), 'Trial');
        load(fullfile(trial_dir, "cm.mat"))
        
        if ~isfile(fullfile(trial_dir, "mu_" + channel + ".mat")), continue, end
        load(fullfile(trial_dir, "mu_" + channel + ".mat"))
        
        cm = cm(1:size(ua, 1));
        
        myplot(t, ua, cm)
        title(selectivity + " Selectivity PSTH")
        subtitle(upper(channel) + ", Folder: " + datestr(folder_name) + ...
            ", Unit: MU")
        touch(gca)
        clear ua
        
        saveas(gcf, fullfile(figure_dir, lower(channel) + num2str(isession) + "-mu.png"))
        close gcf
        
        singleunits = ls(fullfile(trial_dir, "su_" + channel + "*" + ".mat"));
        nsingleunit = size(singleunits, 1);
        for isingleunit = 1:nsingleunit
            load(fullfile(trial_dir, singleunits(isingleunit, :)))
            myplot(t, ua, cm)
            title(selectivity + " Selectivity PSTH")
            subtitle(upper(channel) + ", Folder: " + datestr(folder_name) + ...
                ", Unit: SU-" + num2str(isingleunit))
            touch(gca)
            clear ua
            
            saveas(gcf, fullfile(figure_dir, lower(channel) + num2str(isession) + ...
                "-su-" + num2str(isingleunit) + ".png"))
            close gcf
        end
        clear singleunits nsingleunit isingleunit
    end
end

%%
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