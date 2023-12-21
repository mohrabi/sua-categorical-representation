clc
clear

addpath('tdms')

overwrite = true;

f = fopen('G:\Codes\Preprocessing\Log\tdms2mat_log-eeg.txt', 'at+');
fprintf(f, strcat(datestr(datetime('now')), '\n'));
fprintf(f, "\tConverting tdms files into mat format...\n");
fprintf(f, "\tOverwrite: " + num2str(overwrite) + "\n");

tic
% for monkey = ["Jenab", "Zebel"]
% monkey_dir = fullfile('G:\Data\Slow', monkey);
monkey = "Zebel";
monkey_dir = "J:\Zebel\Fast";
session_dates = ls(fullfile(monkey_dir, '20*'));
%%
for i_session = 1:size(session_dates, 1)
    clc, disp("Converting tdms filed of monkey " + monkey + ...
        ", session " + num2str(i_session) + "/" + num2str(size(session_dates, 1)))
    fprintf(f, "\tConverting tdms filed of monkey " + monkey + ...
        ", session " + num2str(i_session) + "/" + ...
        num2str(size(session_dates, 1)) + "\n");

    tdms_dir = fullfile(monkey_dir, session_dates(i_session, :), 'tdms');
    mat_dir = fullfile(monkey_dir, session_dates(i_session, :), 'mat');

    if isfolder(mat_dir) && overwrite, rmdir(mat_dir, 's'), end
    mkdir(mat_dir)

    try
        parts = ls(fullfile(tdms_dir, 'channelik4*.tdms'));
        spk_it = [];
        for part = 1:size(parts, 1)
            matFileName = simpleConvertTDMS(fullfile(tdms_dir, parts(part, :)));
            load(matFileName{1})
            spk_it = [spk_it; channelik4Spike_IT.Data];
            delete(matFileName{1})
        end
        save(fullfile(tdms_dir, "..", "mat", "spk_it.mat"), 'spk_it', '-v7.3')
    catch ME
        fprintf(f, "\t\tError: " + ME.message + "\n");
    end

    try
        parts = ls(fullfile(tdms_dir, 'channelik5*.tdms'));
        lfp_it = [];
        for part = 1:size(parts, 1)
            matFileName = simpleConvertTDMS(fullfile(tdms_dir, parts(part, :)));
            load(matFileName{1})
            lfp_it = [lfp_it; channelik5LFP_IT.Data];
            delete(matFileName{1})
        end
        save(fullfile(tdms_dir, "..", "mat", "lfp_it.mat"), 'lfp_it', '-v7.3')
    catch ME
        fprintf(f, "\t\tError: " + ME.message + "\n");
    end

    try
        parts = ls(fullfile(tdms_dir, 'channelik8*.tdms'));
        spk_pfc = [];
        for part = 1:size(parts, 1)
            matFileName = simpleConvertTDMS(fullfile(tdms_dir, parts(part, :)));
            load(matFileName{1})
            spk_pfc = [spk_pfc; channelik8Spike_PFC.Data];
            delete(matFileName{1})
        end
        save(fullfile(tdms_dir, "..", "mat", "spk_pfc.mat"), 'spk_pfc', '-v7.3')
    catch ME
        fprintf(f, "\t\tError: " + ME.message + "\n");
    end

    try
        parts = ls(fullfile(tdms_dir, 'channelik9*.tdms'));
        lfp_pfc = [];
        for part = 1:size(parts, 1)
            matFileName = simpleConvertTDMS(fullfile(tdms_dir, parts(part, :)));
            load(matFileName{1})
            lfp_pfc = [lfp_pfc; channelik9LFP_PFC.Data];
            delete(matFileName{1})
        end
        save(fullfile(tdms_dir, "..", "mat", "lfp_pfc.mat"), 'lfp_pfc', '-v7.3')
    catch ME
        fprintf(f, "\t\tError: " + ME.message + "\n");
    end

    try
        parts = ls(fullfile(tdms_dir, 'channelik12*.tdms'));
        eye_x = [];
        for part = 1:size(parts, 1)
            matFileName = simpleConvertTDMS(fullfile(tdms_dir, parts(part, :)));
            load(matFileName{1})
            eye_x = [eye_x; channelik12ChA.Data];
            delete(matFileName{1})
        end

        parts = ls(fullfile(tdms_dir, 'channelik13*.tdms'));
        eye_y = [];
        for part = 1:size(parts, 1)
            matFileName = simpleConvertTDMS(fullfile(tdms_dir, parts(part, :)));
            load(matFileName{1})
            eye_y = [eye_y; channelik13ChB.Data];
            delete(matFileName{1})
        end

        parts = ls(fullfile(tdms_dir, 'channelik14*.tdms'));
        eye_p = [];
        for part = 1:size(parts, 1)
            matFileName = simpleConvertTDMS(fullfile(tdms_dir, parts(part, :)));
            load(matFileName{1})
            eye_p = [eye_p; channelik14ChC.Data];
            delete(matFileName{1})
        end

        eye_xyp = [eye_x, eye_y, eye_p];
        save(fullfile(tdms_dir, "..", "mat", "eye_xyp.mat"), 'eye_xyp', '-v7.3')
    catch ME
        fprintf(f, "\t\tError: " + ME.message + "\n");
    end

    try
        parts = ls(fullfile(tdms_dir, 'eventik_*.tdms'));
        evt = [];
        tim = [];
        for part = 1:size(parts, 1)
            matFileName = simpleConvertTDMS(fullfile(tdms_dir, parts(part, :)));
            load(matFileName{1})
            evt = [evt; eventikSamples.Data];
            tim = [tim; eventikEvents.Data];
            delete(matFileName{1})
        end
        save(fullfile(tdms_dir, "..", "mat", "evt.mat"), 'evt', 'tim', '-v7.3')
    catch ME
        fprintf(f, "\t\tError: " + ME.message + "\n");
    end
end
% end

fprintf(f, "Finished successfully after " + num2str(round(toc)) + ...
    " seconds.\n\n");
fclose(f);