close all
clear
clc

data_root = 'G:\Data\Fast\Jenab';
sessions = string(ls(fullfile(data_root, '20*')));
nsession = length(sessions);

%%
tic
for isession = 72:nsession
    % Make Directories
    out_dir = fullfile('G:\Codes\Processing\out', 'phs', sessions(isession));
    mkdir(out_dir)
    
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
    
    % Extract Phase
    nfreq = 70;
    nsamp = size(it, 2);
    cfs_it = complex(NaN(ntrial, nfreq, nsamp));
    cfs_pf = complex(NaN(ntrial, nfreq, nsamp));
    
    fb = cwtfilterbank('SignalLength', size(it, 2), ...
        'SamplingFrequency', 1000, ...
        'FrequencyLimits', [1, 120]);
    for itrial = 1:ntrial
        [cfs_it(itrial, :, :), ~] = wt(fb, it(itrial, :));
        [cfs_pf(itrial, :, :), freqs] = wt(fb, pf(itrial, :));
    end
    ph_it = angle(cfs_it(:, :, 2000:end-2000-1));
    am_it = abs(cfs_it(:, :, 2000:end-2000-1));
    ph_pf = angle(cfs_pf(:, :, 2000:end-2000-1));
    am_pf = abs(cfs_pf(:, :, 2000:end-2000-1));
    
    % Plot Phase Distribution On Time-Frequqency Plane
    figure('Units', 'centimeters', 'Position', [0, 0, 20, 36])
    tl = tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');
    colormap('bwr')
    nexttile
    contourf(-200:699, freqs, squeeze(mean(ph_it, 1)), 400, 'LineColor', 'none')
    title("Mean"), axis('tight'), colorbar, set(gca, 'YScale', 'log')
    nexttile
    contourf(-200:699, freqs, squeeze(std(ph_it, 1)), 400, 'LineColor', 'none')
    title("STD"), axis('tight'), colorbar, set(gca, 'YScale', 'log')
    title(tl, "Phase of Local Field Potentials, Aligned to Stimulus Onset")
    saveas(gcf, fullfile(out_dir, 'phase-distribution-tf.png'))
    close gcf
    
    % 
    
%     break
end
toc
%%
% winStartIndices = 1:20:801;
% pIaPc = NaN([size(ph_it), length(winStartIndices)]);
% pPaIc = NaN([size(ph_it), length(winStartIndices)]);
% pIaIc = NaN([size(ph_it), length(winStartIndices)]);
% pPaPc = NaN([size(ph_it), length(winStartIndices)]);
% 
% n = 100;
% 
% tic
% for itrial = 1:ntrial
%     for jphas = 1:length(freqs)
%         for kampl = 1:length(freqs)
%             for ttime = 1:length(winStartIndices)
%                 tind = winStartIndices(ttime);
%                 
%                 aI = squeeze(am_it(itrial, kampl, tind:(tind+n-1)));
%                 pI = squeeze(ph_it(itrial, jphas, tind:(tind+n-1)));
%                 aP = squeeze(am_pf(itrial, kampl, tind:(tind+n-1)));
%                 pP = squeeze(ph_pf(itrial, jphas, tind:(tind+n-1)));
%                 
%                 pIaPc(ktrial, jAmp, iPhase, ttime) = abs(1/n * ...
%                     sum(aP * exp(1j * pI')));
%                 pPaIc(ktrial, jAmp, iPhase, ttime) = abs(1/n * ...
%                     sum(aI * exp(1j * pP')));
%                 pIaIc(ktrial, jAmp, iPhase, ttime) = abs(1/n * ...
%                     sum(aI * exp(1j * pI')));
%                 pPaPc(ktrial, jAmp, iPhase, ttime) = abs(1/n * ...
%                     sum(aP * exp(1j * pP')));
% 
% %                 parfor tPerm = 1:nPerm
% %                     x2 = circshift(x, randi(200));
% %                     pac0(jAmp, iPhase, kTrial, tPerm) = abs(1/n * ...
% %                         sum(x2 * exp(1j * y')));
% %                 end
%                 break
%             end
%             break
%         end
%         break
%     end
%     break
% end
% toc
% 
% %%
% nPerm = 100;
% pac = nan(nBand, 11, size(fp, 1));
% pac0 = nan(nBand, 11, size(fp, 1), nPerm);
% for iPhase = 1:11
%     for jAmp = 1:nBand
%         for kTrial = 1:size(fp, 1)
%             n = length(4250:4750);
%             x = squeeze(abs(hlbrt(kTrial, 4250:4750, jAmp)));
%             y = squeeze(angle(hlbrt(kTrial, 4250:4750, iPhase)));
%             pac(jAmp, iPhase, kTrial) = abs(1/n * ...
%                 sum(x * exp(1j * y')));
%             
%             parfor tPerm = 1:nPerm
%                 x2 = circshift(x, randi(200));
%                 pac0(jAmp, iPhase, kTrial, tPerm) = abs(1/n * ...
%                     sum(x2 * exp(1j * y')));
%             end
%         end
%     end
% end