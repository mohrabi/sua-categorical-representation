close all
clear
clc

monkey = "Zebel";
fast_n = string(ls(fullfile("G:\Data\Fast", monkey, '20*')));
corr_n = string(ls(fullfile("J:\Data\Fast\Zebel\Corrupted", '20*')));
slow_n = string(ls(fullfile("G:\Data\Slow", monkey, '20*')));

fast = datetime(fast_n, 'InputFormat', "yyyy-MM-dd_HH-mm");
slow = datetime(slow_n, 'InputFormat', "yyyy-MM-dd_HH-mm");
corr = datetime(corr_n, 'InputFormat', "yyyy-MM-dd_HH-mm");

sel = NaN(length(slow), 1);
dif = [];
for islow = 1:length(slow)
    x = fast - slow(islow);
    arg = find(x < 0);
    clear x
    
    y = corr - slow(islow);
    argc = find(y < 0);
    clear y
    
    if arg(end) >= length(fast)
        continue
    end
    if ~isempty(argc) && (slow(islow) - corr(argc(end))) < (slow(islow) - fast(arg(end)))
        disp(islow)
    end
    
    sel(islow) = arg(end);
    dif = [dif; slow(islow) - fast(arg(end))];
end
clear slow_no arg argc slow fast corr islow

pair = ["Slow", "Fast"];
for i = 1:length(slow_n)
    if isnan(sel(i)), continue, end
    pair = [pair; slow_n(i), fast_n(sel(i))];
end
clear i sel

writematrix(pair, lower(monkey) + "-paired.csv")

clear monkey dif