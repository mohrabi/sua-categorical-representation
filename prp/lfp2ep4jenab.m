clear
src_root = 'J:\Data\Fast\Jenab';
dst_root = 'G:\Data\Fast\Jenab';

sessions = string(ls(fullfile(src_root, '20*')));
nsession = length(sessions);
for isession = 1:nsession
    src = fullfile(src_root, sessions(isession), 'lfp');
    dst = fullfile(dst_root, sessions(isession), 'Trial');
    
    copyfile(fullfile(src, 'it.mat'), fullfile(dst, 'l_it.mat'));
    copyfile(fullfile(src, 'pfc.mat'), fullfile(dst, 'l_pfc.mat'));
end