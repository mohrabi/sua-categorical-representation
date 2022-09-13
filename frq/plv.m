% Phase Locking Analysis
close all
clear
clc

fb = cwtfilterbank('SignalLength', 900, ...
    'SamplingFrequency', 1000, ...
    'FrequencyLimits', [1, 120]);

% Jenab
sessions = string(ls("G:\Data\Fast\Jenab\20*"));
nsession = length(sessions); 

for isession = 1:nsession
    disp([isession, nsession])
    itc = load(fullfile("G:\Data\Fast\Jenab", sessions(isession), ...
        'Trial', 'l_it.mat')).data;
    pfc = load(fullfile("G:\Data\Fast\Jenab", sessions(isession), ...
        'Trial', 'l_pfc.mat')).data;
    itc(any(isnan(itc), 2), :) = [];
    pfc(any(isnan(pfc), 2), :) = [];
    cm = load(fullfile("G:\Data\Fast\Jenab", sessions(isession), ...
        'Trial', 'cm.mat')).cm;
    ntrial = min([size(itc, 1), size(pfc, 1), length(cm)]);
    cm = cm(1:ntrial);
    itc = itc(1:ntrial, :, :);
    pfc = pfc(1:ntrial, :, :);
    itc = itc - mean(itc, 1);
    pfc = pfc - mean(pfc, 1);
    
    [f, lv] = get_plv(itc, pfc, fb);
    out_dir = fullfile("G:\Codes\Processing\out\plv\jenab-fast", sessions(isession));
    if ~isfolder(out_dir), mkdir(out_dir), end
    save(fullfile(out_dir, 'plv.mat'), 'f', 'lv');
    
    clear f lv
%     break
end

%%
t = -150:10:649;
mu = zeros(1, 51, 80);
muc = zeros(4, 137, 51, 80);

for isession = 1:137
    cm = load(fullfile("G:\Data\Fast\Jenab", sessions(isession), ...
        'Trial', 'cm.mat')).cm;
    load(fullfile("G:\Codes\Processing\out\plv\jenab-fast", ...
        sessions(isession), 'plv.mat'))
    bl = mean(lv(:, :, (t < 0) & (t >= -50)), 3);
    lv = lv - bl;
    
    mu = mu + mean(lv, 1) / 138;
    
    lbl = grablabels('face-body');
    ind = ismember(cm, find(lbl == 0));
    muc(1, isession, :, :) = mean(lv(ind(1:size(lv, 1)), :, :), 1);
    ind = ismember(cm, find(lbl == 1));
    muc(2, isession, :, :) = mean(lv(ind(1:size(lv, 1)), :, :), 1);
    lbl = grablabels('artificial-natural');
    ind = ismember(cm, find(lbl == 0));
    muc(3, isession, :, :) = mean(lv(ind(1:size(lv, 1)), :, :), 1);
    ind = ismember(cm, find(lbl == 1));
    muc(4, isession, :, :) = mean(lv(ind(1:size(lv, 1)), :, :), 1);
end

%%
figure('Units', 'centimeters', 'Position', [0, 0, 36, 15])

mucm = squeeze(mean(muc(:, :, (f <= 70) & (f > 50), :), 3));
plot(t, nanmean(squeeze(mucm(4, :, :))))
hold on
mucm = squeeze(mean(muc(:, :, (f <= 32) & (f > 13), :), 3));
plot(t, nanmean(squeeze(mucm(4, :, :))))
% plotmetric(t, {squeeze(mucm(1, :, :))})
legend("Gamma", "Beta")
title("Inanimate")
xlim([-50, 500])
touch(gca)

%%
figure('Units', 'centimeters', 'Position', [0, 0, 36, 15])

% mucm = squeeze(mean(muc(:, :, (f <= 70) & (f > 50), :), 3));
mucm = squeeze(mean(muc(:, :, (f <= 32) & (f > 13), :), 3));
plot(t, nanmean(squeeze(mucm(1, :, :))))
hold on
plot(t, nanmean(squeeze(mucm(4, :, :))))
% plotmetric(t, {squeeze(mucm(1, :, :))})
legend("Face", "Inanimate")
title("Beta")
xlim([-50, 500])
touch(gca)

%%
subplot(2, 2, 1)
tInd = (t > -50) & (t < 500);
fInd = (f > 13) & (f <= 120);
X = squeeze(mean(muc(1, :, fInd, tInd), 2));
contourf(t(tInd), f(fInd), X , 400, 'LineColor', 'none')
set(gca, 'YScale', 'log')
colormap('jet')
colorbar
title('Face')

subplot(2, 2, 2)
X = squeeze(mean(muc(2, :, fInd, tInd), 2));
contourf(t(tInd), f(fInd), X , 400, 'LineColor', 'none')
set(gca, 'YScale', 'log')
colormap('jet')
colorbar
title('Body')

subplot(2, 2, 3)
X = squeeze(nanmean(muc(3, :, fInd, tInd), 2));
contourf(t(tInd), f(fInd), X , 400, 'LineColor', 'none')
set(gca, 'YScale', 'log')
colormap('jet')
colorbar
title('Artificial')

subplot(2, 2, 4)
X = squeeze(nanmean(muc(4, :, fInd, tInd), 2));
contourf(t(tInd), f(fInd), X , 400, 'LineColor', 'none')
set(gca, 'YScale', 'log')
colormap('jet')
colorbar
title('Natural')

%%
figure('Units', 'centimeters', 'Position', [0, 0, 36, 15])

addpath('C:\toolbox\Separation Toolbox\')
mucm = squeeze(mean(muc(:, :, (f <= 70) & (f > 60), :), 3));
plotmetric(t, {squeeze(mucm(1, :, :)), ...
    squeeze(mucm(2, :, :)), ...
    squeeze(mucm(3, :, :)), ...,
    squeeze(mucm(4, :, :))}, 3)
legend("Face", "Body", "Artificial", "Natural")
xlim([-50, 500])
touch(gca)

%%
subplot(1, 2, 1)
tInd = (t > 0) & (t < 500);
fInd = (f > 13) & (f <= 32);
contourf(t(tInd), f(fInd), squeeze(mu(:, fInd, tInd)), 400, 'LineColor', 'none')
set(gca, 'YScale', 'log')
colormap('jet')
colorbar

subplot(1, 2, 2)
tInd = (t > 0) & (t < 500);
fInd = (f > 32) & (f <= 120);
contourf(t(tInd), f(fInd), squeeze(mu(:, fInd, tInd)), 400, 'LineColor', 'none')
set(gca, 'YScale', 'log')
colormap('jet')
colorbar

%%
% mun = (mu - min(mu, [], 'all')) / (max(mu, [], 'all') - min(mu, [], 'all'));
mun = mu;

%%
axs = NaN(3, 1);
fInd = (f >= 13) & (f < 32);

axs(1) = subplot(1, 3, 1);
tInd = (t >= -100) & (t < 0);
[c, cc] = contourf(t(tInd), f(fInd), squeeze(mun(:, fInd, tInd)), 400, 'LineColor', 'none');
set(gca, 'YScale', 'log')
colormap('jet')
colorbar
% caxis([-20e-3, 3e-3])

axs(2) = subplot(1, 3, 2);
tInd = (t >= 100) & (t < 200);
contourf(t(tInd), f(fInd), squeeze(mun(:, fInd, tInd)), 400, 'LineColor', 'none')
set(gca, 'YScale', 'log')
colormap('jet')
colorbar
% caxis([-20e-3, 3e-3])

axs(3) = subplot(1, 3, 3);
tInd = (t >= 200) & (t < 300);
contourf(t(tInd), f(fInd), squeeze(mun(:, fInd, tInd)), 400, 'LineColor', 'none')
set(gca, 'YScale', 'log')
colormap('jet')
colorbar
caxis([-20e-3, 3e-3])

%%
figure

nexttile
lvm = lv(:, (f>13) & (f<=32), :);
contourf(-150:10:649, f((f>13) & (f<=32)), squeeze(mean(lvm, 1)), 100, 'LineColor', 'none')
colormap(jet)
colorbar

nexttile
lvm = lv(:, (f>32) & (f<=128), :);
contourf(-150:10:649, f((f>32) & (f<=128)), squeeze(mean(lvm, 1)), 100, 'LineColor', 'none')
colormap(jet)
colorbar

%%

monkey = "Jenab";
selectivity = "Slow";

addpath(genpath('C:\toolbox\Chronux'))
addpath('C:\toolbox\Separation Toolbox')
addpath('G:\Codes\Common')
addpath('src')

data_root = fullfile('G:\Data', selectivity, monkey);
sessions = string(ls(fullfile(data_root, '20*')));
nsession = length(sessions);

if selectivity == "Fast"
    t = -200:699;
elseif selectivity == "Slow"
    t = -200:999;
end
fs = 1000;

plv1 = []; plv2 = []; plv3 = [];
for isession = 2:nsession
    clc, disp([isession, nsession])
    [it, pf, cm] = load_data(fullfile(data_root, sessions(isession)));
    it = it - mean(it, 1);
    pf = pf - mean(pf, 1);
    
    nfreq = 55;
    [ntrial, nsamp] = size(it);
    cfs_it = complex(NaN(ntrial, nfreq, nsamp));
    cfs_pf = complex(NaN(ntrial, nfreq, nsamp));
    
    fb = cwtfilterbank('SignalLength', size(it, 2), ...
        'SamplingFrequency', 1000, ...
        'FrequencyLimits', [1, 120]);
    for itrial = 1:ntrial
        [cfs_it(itrial, :, :), ~] = wt(fb, it(itrial, :));
        [cfs_pf(itrial, :, :), freqs] = wt(fb, pf(itrial, :));
    end
    
    trange = [70, 270];
    pit = angle(cfs_it(:, :, (t >= trange(1)) & (t < trange(2))));
    ppf = angle(cfs_pf(:, :, (t >= trange(1)) & (t < trange(2))));
    X = abs(mean(exp(1j*(ppf - pit)), 3));

    trange = [270, 470];
    pit = angle(cfs_it(:, :, (t >= trange(1)) & (t < trange(2))));
    ppf = angle(cfs_pf(:, :, (t >= trange(1)) & (t < trange(2))));
    Y = abs(mean(exp(1j*(ppf - pit)), 3));
    
    trange = [-130, 70];
    pit = angle(cfs_it(:, :, (t >= trange(1)) & (t < trange(2))));
    ppf = angle(cfs_pf(:, :, (t >= trange(1)) & (t < trange(2))));
    Z = abs(mean(exp(1j*(ppf - pit)), 3));
    
    plv1 = [plv1; mean(Z, 1)];
    plv2 = [plv2; mean(X, 1)];
    plv3 = [plv3; mean(Y, 1)];
end
clear isession nsession sessions
% save(fullfile('out', 'cca', 'early_vs_late_overall.mat'), 'f1', 'f2', 'Ct1', 'Ct2')

%%
plot(freqs, mean(plv1, 1), 'LineWidth', 3)
hold on
plot(freqs, mean(plv2, 1), 'LineWidth', 3)
plot(freqs, mean(plv3, 1), 'LineWidth', 3)

%%
% load(fullfile('out', 'cca', 'early_vs_late_overall.mat'), 'f1', 'f2', 'Ct1', 'Ct2')
figure('Units', 'centimeters', 'Position', [0, 0, 36, 15])
subplot(1, 4, 1)
plot(freqs, mean(X, 1), 'LineWidth', 3)
hold on
plot(freqs, mean(Y, 1), 'LineWidth', 3)
plot(freqs, mean(Z, 1), 'LineWidth', 3)
xlim([1, 30])
touch(gca)
subplot(1, 4, 2:4)
plot(freqs, mean(X, 1), 'LineWidth', 3)
hold on
plot(freqs, mean(Y, 1), 'LineWidth', 3)
plot(freqs, mean(Z, 1), 'LineWidth', 3)
xlim([30, 120])
touch(gca)


%% 
function [it, pf, cm] = load_data(path)
it = load(fullfile(path, "Trial/l_it.mat")).data;
pf = load(fullfile(path, "Trial/l_pfc.mat")).data;
cm = load(fullfile(path, "Trial/cm.mat")).cm;
ntrial = min([length(cm), size(it, 1), size(pf, 1)]);
it = it(1:ntrial, :);
pf = pf(1:ntrial, :);
cm = cm(1:ntrial);
ind = isnan(cm) | all(isnan(it), 2) | all(isnan(pf), 2);
cm = cm(~ind);
it = it(~ind, :);
pf = pf(~ind, :);
ntrial = length(cm);
end

function [C, P, f] = coherency(X, Y, fs)
assert(all(size(X) == size(Y)));

tapers = dpss(size(X, 2), 3, 5);
tapers = tapers * sqrt(fs);

frange = [2, 120];

[~, nsample] = size(X);
nfft = 2 ^ (nextpow2(nsample));

FX = mtfftc(X', tapers, nfft, 1000);
FY = mtfftc(Y', tapers, nfft, 1000);

SXY = squeeze(mean(conj(FX) .* FY, 2));
SX = squeeze(mean(conj(FX) .* FX, 2));
SY = squeeze(mean(conj(FY) .* FY, 2));

SX = mean(SX, 2); SY = mean(SY, 2); SXY = mean(SXY, 2);

C = SXY ./sqrt(SX .* SY);
f = linspace(0, fs, nfft);

P = angle(C((f < frange(2)) & (f >= frange(1))));
C = abs(C((f < frange(2)) & (f >= frange(1))));
f = f((f < frange(2)) & (f >= frange(1)));
end

function plotmetric(bt, si, lw)
ctbl = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330;
    0.6350 0.0780 0.1840];

hold on
for ii = 1:length(si)
    si_no = si{ii};
    si_no(~any(isfinite(si_no), 2), :) = [];
    data = mean(si_no, 1);
    % conf = prctile(si_no, [2.5, 97.5], 1);
    conf = data + [-1.96; 1.96] * std(si_no, [], 1) / sqrt(size(si_no, 1));
    plot(bt, data, 'LineWidth', 2, 'Color', ctbl(ii, :), 'LineWidth', lw)
    fill([bt, flip(bt)], [conf(1, :), flip(conf(2, :))], ctbl(ii, :), ...
        'EdgeColor', 'none', 'FaceAlpha', .1, 'HandleVisibility', 'off')
end
xlim([bt(1), bt(end)])
end

function [freqs, lv] = get_plv(X1, X2, fb)

[nTrial, nsample] = size(X1);
nfreq = 51;

cfs_itc = complex(NaN(nTrial, nfreq, nsample));
cfs_pfc = complex(NaN(nTrial, nfreq, nsample));
for iTrial = 1:nTrial
    [cfs_itc(iTrial, :, :), ~]      = wt(fb, X1(iTrial, :));
    [cfs_pfc(iTrial, :, :), freqs]  = wt(fb, X2(iTrial, :));
end

w0 = 1:10:800;
% nPerm = 200;

lv  = NaN(nTrial, nfreq, 80);
% lv0 = NaN(nTrial, nfreq, 80, nPerm);
for iTrial = 1:nTrial
%     disp([iTrial, nTrial])
    for jWin = 1:length(w0)
        p1 = squeeze(angle(cfs_itc(iTrial, :, :)));
        p2 = squeeze(angle(cfs_pfc(iTrial, :, :)));
        
        w = w0(jWin) : (w0(jWin) + 99);
        lv(iTrial, :, jWin) = abs(mean(exp(1j * (p1(:, w) - p2(:, w))), 2));
        
%         for kPerm = 1:nPerm
%             p1 = circshift(p1, randi(99), 2);
%             lv0(iTrial, :, jWin, kPerm) = abs(mean(exp(1j * (p1(:, w) - p2(:, w))), 2));
%         end
    end
end

% t = -150:10:649;
% bl = mean(lv(:, :, t < 0), [1, 3]);
% lv = lv - bl;
% contourf(t, freqs, squeeze(mean(lv, 1)), 400, 'LineColor', 'none')
% set(gca, 'YScale', 'log')
% lv = (lv - nanmean(lv0, 4)) ./ nanstd(lv0, [], 4);

end