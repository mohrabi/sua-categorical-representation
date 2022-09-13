% Phase of Cross-Coherency Analysis
close all
clear
clc

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

Pt1 = []; Pt2 = [];
for isession = 1:nsession
    [it, pf, cm] = load_data(fullfile(data_root, sessions(isession)));
    trange = [70, 270];
    [~, P1, f1] = coherency(it(:, (t >= trange(1)) & (t < trange(2))), ...
        pf(:, (t >= trange(1)) & (t < trange(2))), fs);
    trange = [270, 470];
    [~, P2, f2] = coherency(it(:, (t >= trange(1)) & (t < trange(2))), ...
        pf(:, (t >= trange(1)) & (t < trange(2))), fs);
    
    Pt1 = [Pt1; P1'];
    Pt2 = [Pt2; P2'];
end
clear isession nsession sessions
% save(fullfile('out', 'cca', 'early_vs_late_overall.mat'), 'f1', 'f2', 'Ct1', 'Ct2')

%%
plot(f1, angle(mean(exp(1j*Pt1)))')
hold on
plot(f2, angle(mean(exp(1j*Pt2)))')


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
    si_no(:, ~any(isfinite(si_no), 2)) = [];
    data = mean(si_no, 1);
    % conf = prctile(si_no, [2.5, 97.5], 1);
    conf = data + [-1; 1] * std(si_no, [], 1);
    plot(bt, data, 'LineWidth', 2, 'Color', ctbl(ii, :), 'LineWidth', lw)
    fill([bt, flip(bt)], [conf(1, :), flip(conf(2, :))], ctbl(ii, :), ...
        'EdgeColor', 'none', 'FaceAlpha', .1, 'HandleVisibility', 'off')
end
xlim([bt(1), bt(end)])
end
