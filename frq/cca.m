% Cross-Coherence Analysis
close all
clear
clc

addpath(genpath('C:\toolbox\Chronux'))
addpath('C:\toolbox\Separation Toolbox')
addpath('G:\Codes\Common')
addpath('src')

data_root = 'G:\Data\Slow\Jenab';
sessions = string(ls(fullfile(data_root, '20*')));
nsession = length(sessions);

t = -200:999;
fs = 1000;

Ct1 = []; Ct2 = [];
for isession = 1:nsession
    [it, pf, cm] = load_data(fullfile(data_root, sessions(isession)));
    trange = [70, 270];
    [C1, f1] = coherency(it(:, (t >= trange(1)) & (t < trange(2))), ...
        pf(:, (t >= trange(1)) & (t < trange(2))), fs);
    trange = [270, 470];
    [C2, f2] = coherency(it(:, (t >= trange(1)) & (t < trange(2))), ...
        pf(:, (t >= trange(1)) & (t < trange(2))), fs);
    
    Ct1 = [Ct1; C1'];
    Ct2 = [Ct2; C2'];
end
clear isession nsession sessions

% save(fullfile('out', 'cca', 'early_vs_late_overall.mat'), 'f1', 'f2', 'Ct1', 'Ct2')
%%
clear
clc
load(fullfile('out', 'cca', 'early_vs_late_overall.mat'), 'f1', 'f2', 'Ct1', 'Ct2')
figure('Units', 'centimeters', 'Position', [0, 0, 36, 15])
subplot(1, 4, 1)
plot(f1, mean(Ct1, 1), 'LineWidth', 3)
hold on
plot(f2, mean(Ct2, 1), 'LineWidth', 3)
xlim([1, 30])
touch(gca)

subplot(1, 4, 2:4)
plot(f1, mean(Ct1, 1), 'LineWidth', 3)
hold on
plot(f2, mean(Ct2, 1), 'LineWidth', 3)
xlim([30, 120])
touch(gca)
%%
close all
clear
clc

addpath(genpath('C:\toolbox\Chronux'))
addpath('C:\toolbox\Separation Toolbox')
addpath('G:\Codes\Common')
addpath('src')

data_root = 'G:\Data\Fast\Jenab';
sessions = string(ls(fullfile(data_root, '20*')));
nsession = length(sessions);

t = -200:999;
fs = 1000;

Ct1.face = []; Ct1.body = []; Ct1.animate = []; Ct1.inanimate = [];
Ct2.face = []; Ct2.body = []; Ct2.animate = []; Ct2.inanimate = [];
for isession = 1:nsession
    disp([isession, nsession])
    [it, pf, cm] = load_data(fullfile(data_root, sessions(isession)));
    trange = [70, 270];
    lbl = grablabels('super-ordinate');
    cond = find(lbl==0);
    [C1, ~] = coherency(it(ismember(cm, cond), (t >= trange(1)) & (t < trange(2))), ...
        pf(ismember(cm, cond), (t >= trange(1)) & (t < trange(2))), fs);
    Ct1.animate = [Ct1.animate; C1'];
    cond = find(lbl==1);
    [C1, ~] = coherency(it(ismember(cm, cond), (t >= trange(1)) & (t < trange(2))), ...
        pf(ismember(cm, cond), (t >= trange(1)) & (t < trange(2))), fs);
    Ct1.inanimate = [Ct1.inanimate; C1'];
    lbl = grablabels('face-body');
    cond = find(lbl==0);
    [C1, ~] = coherency(it(ismember(cm, cond), (t >= trange(1)) & (t < trange(2))), ...
        pf(ismember(cm, cond), (t >= trange(1)) & (t < trange(2))), fs);
    Ct1.face = [Ct1.face; C1'];
    cond = find(lbl==1);
    [C1, f1] = coherency(it(ismember(cm, cond), (t >= trange(1)) & (t < trange(2))), ...
        pf(ismember(cm, cond), (t >= trange(1)) & (t < trange(2))), fs);
    Ct1.body = [Ct1.body; C1'];
    
    trange = [270, 470];
    lbl = grablabels('super-ordinate');
    cond = find(lbl==0);
    [C2, ~] = coherency(it(ismember(cm, cond), (t >= trange(1)) & (t < trange(2))), ...
        pf(ismember(cm, cond), (t >= trange(1)) & (t < trange(2))), fs);
    Ct2.animate = [Ct2.animate; C2'];
    cond = find(lbl==1);
    [C2, ~] = coherency(it(ismember(cm, cond), (t >= trange(1)) & (t < trange(2))), ...
        pf(ismember(cm, cond), (t >= trange(1)) & (t < trange(2))), fs);
    Ct2.inanimate = [Ct2.inanimate; C2'];
    lbl = grablabels('face-body');
    cond = find(lbl==0);
    [C2, ~] = coherency(it(ismember(cm, cond), (t >= trange(1)) & (t < trange(2))), ...
        pf(ismember(cm, cond), (t >= trange(1)) & (t < trange(2))), fs);
    Ct2.face = [Ct2.face; C2'];
    cond = find(lbl==1);
    [C2, f2] = coherency(it(ismember(cm, cond), (t >= trange(1)) & (t < trange(2))), ...
        pf(ismember(cm, cond), (t >= trange(1)) & (t < trange(2))), fs);
    Ct2.body = [Ct2.body; C2'];
end
clear isession nsession sessions

save(fullfile('out', 'cca', 'jenab-slow-early_vs_late_categorical.mat'), 'f1', 'f2', 'Ct1', 'Ct2')


%%
clc
clear
selectivity = "slow";

load(fullfile('out', 'cca', "jenab" + "-" + selectivity + "-early_vs_late_categorical.mat"), ...
    'f1', 'f2', 'Ct1', 'Ct2')

figure('Units', 'centimeters', 'Position', [0, 0, 36, 15])
subplot(1, 4, 1)
plot(f1, mean(Ct1.face, 1), 'LineWidth', 3)
hold on
plot(f2, mean(Ct2.face, 1), 'LineWidth', 3)
xlim([1, 30])
touch(gca)
subplot(1, 4, 2:4)
plot(f1, mean(Ct1.face, 1), 'LineWidth', 3)
hold on
plot(f2, mean(Ct2.face, 1), 'LineWidth', 3)
xlim([30, 120])
legend("Early", "Late")
touch(gca)
suptitle("Face")
saveas(gcf, fullfile('out', 'cca', 'jenab-fast-face.jpg'))
close gcf

figure('Units', 'centimeters', 'Position', [0, 0, 36, 15])
subplot(1, 4, 1)
plot(f1, mean(Ct1.body, 1), 'LineWidth', 3)
hold on
plot(f2, mean(Ct2.body, 1), 'LineWidth', 3)
xlim([1, 30])
touch(gca)
subplot(1, 4, 2:4)
plot(f1, mean(Ct1.body, 1), 'LineWidth', 3)
hold on
plot(f2, mean(Ct2.body, 1), 'LineWidth', 3)
xlim([30, 120])
legend("Early", "Late")
touch(gca)
suptitle("Body")
saveas(gcf, fullfile('out', 'cca', 'jenab-fast-body.jpg'))
close gcf


figure('Units', 'centimeters', 'Position', [0, 0, 36, 15])
subplot(1, 4, 1)
plot(f1, mean(Ct1.animate, 1), 'LineWidth', 3)
hold on
plot(f2, mean(Ct2.animate, 1), 'LineWidth', 3)
xlim([1, 30])
touch(gca)
subplot(1, 4, 2:4)
plot(f1, mean(Ct1.animate, 1), 'LineWidth', 3)
hold on
plot(f2, mean(Ct2.animate, 1), 'LineWidth', 3)
xlim([30, 120])
legend("Early", "Late")
touch(gca)
suptitle("Animate")
saveas(gcf, fullfile('out', 'cca', 'jenab-fast-animate.jpg'))
close gcf

figure('Units', 'centimeters', 'Position', [0, 0, 36, 15])
subplot(1, 4, 1)
plot(f1, mean(Ct1.inanimate, 1, 'omitnan'), 'LineWidth', 3)
hold on
plot(f2, mean(Ct2.inanimate, 1, 'omitnan'), 'LineWidth', 3)
xlim([1, 30])
touch(gca)
subplot(1, 4, 2:4)
plot(f1, mean(Ct1.inanimate, 1, 'omitnan'), 'LineWidth', 3)
hold on
plot(f2, mean(Ct2.inanimate, 1, 'omitnan'), 'LineWidth', 3)
xlim([30, 120])
legend("Early", "Late")
touch(gca)
suptitle("Inanimate")
saveas(gcf, fullfile('out', 'cca', 'jenab-fast-inanimate.jpg'))
close gcf


figure('Units', 'centimeters', 'Position', [0, 0, 36, 15])
subplot(1, 4, 1)
plot(f1, mean(Ct1.inanimate, 1, 'omitnan'), 'LineWidth', 3)
hold on
plot(f1, mean(Ct1.animate, 1, 'omitnan'), 'LineWidth', 3)
plot(f1, mean(Ct1.face, 1, 'omitnan'), 'LineWidth', 3)
plot(f1, mean(Ct1.body, 1, 'omitnan'), 'LineWidth', 3)
xlim([1, 30])
touch(gca)
subplot(1, 4, 2:4)
plot(f1, mean(Ct1.inanimate, 1, 'omitnan'), 'LineWidth', 3)
hold on
plot(f1, mean(Ct1.animate, 1, 'omitnan'), 'LineWidth', 3)
plot(f1, mean(Ct1.face, 1, 'omitnan'), 'LineWidth', 3)
plot(f1, mean(Ct1.body, 1, 'omitnan'), 'LineWidth', 3)
xlim([30, 120])
legend("Inanimate", "Animate", "Face", "Body")
touch(gca)
suptitle("Early")
saveas(gcf, fullfile('out', 'cca', 'jenab-fast-early.jpg'))
close gcf

figure('Units', 'centimeters', 'Position', [0, 0, 36, 15])
subplot(1, 4, 1)
plot(f1, mean(Ct2.inanimate, 1, 'omitnan'), 'LineWidth', 3)
hold on
plot(f1, mean(Ct2.animate, 1, 'omitnan'), 'LineWidth', 3)
plot(f1, mean(Ct2.face, 1, 'omitnan'), 'LineWidth', 3)
plot(f1, mean(Ct2.body, 1, 'omitnan'), 'LineWidth', 3)
xlim([1, 30])
touch(gca)
subplot(1, 4, 2:4)
plot(f1, mean(Ct2.inanimate, 1, 'omitnan'), 'LineWidth', 3)
hold on
plot(f1, mean(Ct2.animate, 1, 'omitnan'), 'LineWidth', 3)
plot(f1, mean(Ct2.face, 1, 'omitnan'), 'LineWidth', 3)
plot(f1, mean(Ct2.body, 1, 'omitnan'), 'LineWidth', 3)
xlim([30, 120])
legend("Inanimate", "Animate", "Face", "Body")
touch(gca)
suptitle("Late")
saveas(gcf, fullfile('out', 'cca', 'jenab-fast-late.jpg'))
close gcf

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

function [C, f] = coherency(X, Y, fs)
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
