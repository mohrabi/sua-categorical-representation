close all
clear
clc

%% Zebel
load("G:\Data\Fast\Zebel\itcNeuralInfo.mat")
load("G:\Data\Fast\Zebel\pfcNeuralInfo.mat")

itc = NaN * zeros(170, size(itcNeuralInfo, 1)-1, 900);
for iRow = 2:size(itcNeuralInfo, 1)
    ua = load(itcNeuralInfo(iRow, 1)).ua;
    cm = load(fullfile(fileparts(itcNeuralInfo(iRow, 1)), 'cm.mat')).cm;
    assert(~any(isnan(ua), 'all'));
    
    for iStm = unique(cm)'
        itc(iStm, iRow-1, :) = mean(ua(cm == iStm, :), 1);
    end
end

pfc = NaN * zeros(170, size(pfcNeuralInfo, 1)-1, 900);
for iRow = 2:size(pfcNeuralInfo, 1)
    ua = load(pfcNeuralInfo(iRow, 1)).ua;
    cm = load(fullfile(fileparts(pfcNeuralInfo(iRow, 1)), 'cm.mat')).cm;
    assert(~any(isnan(ua), 'all'));
    
    for iStm = unique(cm)'
        pfc(iStm, iRow-1, :) = mean(ua(cm == iRow, :), 1);
    end
end

save('G:\Data\Fast\Zebel\ITC.mat', 'itc')
save('G:\Data\Fast\Zebel\PFC.mat', 'pfc')

%% Jenab
load("G:\Data\Fast\Jenab\itcNeuralInfo.mat")
load("G:\Data\Fast\Jenab\pfcNeuralInfo.mat")

itc = NaN * zeros(165, size(itcNeuralInfo, 1)-1, 900);
for iRow = 2:size(itcNeuralInfo, 1)
    ua = load(itcNeuralInfo(iRow, 1)).ua;
    cm = load(fullfile(fileparts(itcNeuralInfo(iRow, 1)), 'cm.mat')).cm;
    assert(~any(isnan(ua), 'all'));
    
    for iStm = unique(cm)'
        itc(iStm, iRow-1, :) = mean(ua(cm == iStm, :), 1);
    end
end

pfc = NaN * zeros(165, size(pfcNeuralInfo, 1)-1, 900);
for iRow = 2:size(pfcNeuralInfo, 1)
    ua = load(pfcNeuralInfo(iRow, 1)).ua;
    cm = load(fullfile(fileparts(pfcNeuralInfo(iRow, 1)), 'cm.mat')).cm;
    assert(~any(isnan(ua), 'all'));
    
    for iStm = unique(cm)'
        pfc(iStm, iRow-1, :) = mean(ua(cm == iRow, :), 1);
    end
end

save('G:\Data\Fast\Jenab\ITC.mat', 'itc')
save('G:\Data\Fast\Jenab\PFC.mat', 'pfc')

%% Both
clear
close all
clc

itc_j = load('G:\Data\Fast\Jenab\ITC.mat', 'itc').itc;
itc_z = load('G:\Data\Fast\Zebel\ITC.mat', 'itc').itc;
pfc_j = load('G:\Data\Fast\Jenab\PFC.mat', 'pfc').pfc;
pfc_z = load('G:\Data\Fast\Zebel\PFC.mat', 'pfc').pfc;

itc = cat(2, itc_j, itc_z(1:165, :, :));
pfc = cat(2, pfc_j, pfc_z(1:165, :, :));

save("G:\Data\Fast\Both\ITC.mat", 'itc')
save("G:\Data\Fast\Both\PFC.mat", 'pfc')
