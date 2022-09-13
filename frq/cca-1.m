% Cross-Coherence Analysis
close all
clear
clc

addpath(genpath('C:\toolbox\Chronux'))
addpath('C:\toolbox\Separation Toolbox')
addpath('G:\Codes\Common')


data_root = 'G:\Data\Fast\Jenab';
sessions = string(ls(fullfile(data_root, '20*')));
nsession = length(sessions);

%%
tic
params.tapers   = [3, 5];
params.pad      = 0;
params.Fs       = 1000;
params.fpass    = [4 120];
params.trialave = 0;
params.err      = [2, .95];

for isession = 1:nsession
    clc, disp([isession, nsession])
    % Make Directories
    out_path = fullfile('G:\Codes\Processing\out', 'cca', sessions(isession));
    if ~isfolder(out_path), mkdir(out_path), end
    
    % Load Session Data
    it = load(fullfile(data_root, sessions(isession), "Trial/l_it.mat")).data;
    pf = load(fullfile(data_root, sessions(isession), "Trial/l_pfc.mat")).data;
    cm = load(fullfile(data_root, sessions(isession), "Trial/cm.mat")).cm;
    ntrial = min([length(cm), size(it, 1), size(pf, 1)]);
    it = it(1:ntrial, :);
    pf = pf(1:ntrial, :);
    cm = cm(1:ntrial);
    ind = isnan(cm) | all(isnan(it), 2) | all(isnan(pf), 2);
    cm = cm(~ind);
    it = it(~ind, :);
    pf = pf(~ind, :);
    ntrial = min([length(cm), size(it, 1), size(pf, 1)]);
    clear ind
    %%
    t = -2200:2699;
    X = it(:, )
    %%
    [C, phi, S12, S1, S2, t, f, confC, phistd, Cerr] = ...
        cohgramc(it', pf', [.3, .025], params);
    t = 1000 * t - 2001;
    %%
    break
    save(fullfile(out_path, 'cc.mat'), 't', 'f', 'cm', 'C', 'Cerr')
    clear out_path
end
clear isession
toc
%%
ifreq = 2;
cc = NaN(nsession, 4, 29, 189);

for isession = 1:nsession
    clc, isession
    load(fullfile('G:\Codes\Processing\out\cca', sessions(isession), "cc.mat"))
    C = C - mean(C((t < 50) & (t >= -100), :, :), 1);
    lbl = grablabels('super-ordinate');
    cond = ismember(cm, find(lbl==0));
    cc(isession, 1, :, :) = mean(C(:,:,cond), 3)';
    cond = ismember(cm, find(lbl==1));
    cc(isession, 2, :, :) = mean(C(:,:,cond), 3)';
    lbl = grablabels('face-body');
    cond = ismember(cm, find(lbl==0));
    cc(isession, 3, :, :) = mean(C(:,:,cond), 3)';
    cond = ismember(cm, find(lbl==1));
    cc(isession, 4, :, :) = mean(C(:,:,cond), 3)';
end 
%%
nanind = any(isnan(cc), [2, 3, 4]);
infind = any(isinf(cc), [2, 3, 4]);
ind = ~nanind & ~infind;
for ifreq = 1:29
    figure
    for i = 1:4
        m = squeeze(mean(cc(ind, i, ifreq, :), 1));
        s = 1.96 * squeeze(std(cc(ind, i, ifreq, :), [], 1)) / sqrt(sum(ind));
        plotmetricx(t(:)', m(:)', [m(:)-s(:), m(:)+s(:)]')
    end
    xlim([-200, 700])
end
%%
X = squeeze(cc(ind, 3, 1, :)) - squeeze(mean(cc(ind, 3, 1, (t<0)&(t>=-200)), 4));
plotmetricx(t, mean(X, 1), prctile(X, [2.5, 97.5]))
hold on
X = squeeze(cc(ind, 4, 1, :)) - squeeze(mean(cc(ind, 4, 1, (t<0)&(t>=-200)), 4));
plotmetricx(t, mean(X, 1), prctile(X, [2.5, 97.5]))

X = squeeze(cc(ind, 1, 1, :)) - squeeze(mean(cc(ind, 1, 1, (t<0)&(t>=-200)), 4));
plotmetricx(t, mean(X, 1), prctile(X, [2.5, 97.5]))
X = squeeze(cc(ind, 2, 1, :)) - squeeze(mean(cc(ind, 2, 1, (t<0)&(t>=-200)), 4));
plotmetricx(t, mean(X, 1), prctile(X, [2.5, 97.5]))
xlim([-200, 700])
%%
for ifreq = 1%:5%length(f)
    figure('Units', 'centimeters', 'Position', [0, 0, 36, 15])
    
    lbl = grablabels('super-ordinate');
    cond = ismember(cm, find(lbl==0));
    plotmetricx(t, squeeze(mean(C(:,ifreq,cond), 3)), squeeze(mean(Cerr(:,:,ifreq,cond), 4)))
    cond = ismember(cm, find(lbl==1));
    plotmetricx(t, squeeze(mean(C(:,ifreq,cond), 3)), squeeze(mean(Cerr(:,:,ifreq,cond), 4)))
    lbl = grablabels('face-body');
    cond = ismember(cm, find(lbl==0));
    plotmetricx(t, squeeze(mean(C(:,ifreq,cond), 3)), squeeze(mean(Cerr(:,:,ifreq,cond), 4)))
    cond = ismember(cm, find(lbl==1));
    plotmetricx(t, squeeze(mean(C(:,ifreq,cond), 3)), squeeze(mean(Cerr(:,:,ifreq,cond), 4)))
    
    xlim([-200, 700])
    title("Cross       -Coherence: Frequency = " + num2str(f(ifreq)))
    legend("Animate", "Inanimate", "Face", "Body")
    touch(gca)
end

function plotmetricx(t, m, c)
ctbl = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330;
    0.6350 0.0780 0.1840];

k = normpdf(-2:2, 0, 2);
k = k / sum(k);
m = conv(m, k, 'same');
c = conv2(c, k, 'same');

ax = gca;
ind = ax.ColorOrderIndex;
hold on
plot(t, m, 'LineWidth', 2, 'Color', ctbl(ind, :), 'LineWidth', 3)
fill([t, flip(t)], [c(1, :), flip(c(2, :))], ctbl(ind, :), ...
    'EdgeColor', 'none', 'FaceAlpha', .1, 'HandleVisibility', 'off')

ind = ind + 1;
if ind > size(ax.ColorOrder, 1)
    ind = 1;
end
ax.ColorOrderIndex = ind;
end
