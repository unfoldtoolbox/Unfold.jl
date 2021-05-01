%% Script to implement new test cases for the Unfold.jl toolbox

% Parameters for the simulation:

% Overlap: -1 without overlap, 0.4 with overlap
% Noise: 0 (maybe changed later)
% Sampling rate (srate): 20 (1000 for long data)
% Data length: 60 minutes (i.e. 60*60) for long data, 10 minutes (default;
% i.e. 10*60) else
% Basis (i.e. kernel shape): 'box'
% Effect sizes:
%   2 for intercept
%   3 for categorical variables
%   4 for continuous variables (range = [-1,1])
%   2 and 3 for intecepts of two events


% set random seed
rng(42);

%% Case 1a:  eventA: y~1 (without overlap)
intercept = struct();
intercept.eventname = 'eventA';
intercept.type = 'intercept';
intercept.predictorName = 'intercept';
intercept.overlap = -1;
intercept.effectsize = 2;

signals{1} = intercept;
%signals{1}.range = nan;

EEG_sim_test_case_1a = simulate_data(signals,'noise',0,'srate',20);
save_data(EEG_sim_test_case_1a,'test_case_1a')

%% Case 1b:  eventA: y~1 (with overlap)
intercept = struct();
intercept.eventname = 'eventA';
intercept.type = 'intercept';
intercept.predictorName = 'intercept';
intercept.overlap = 0.4;
intercept.effectsize = 2;

signals{1} = intercept;
%signals{1}.range = nan;

EEG_sim_test_case_1b = simulate_data(signals,'noise',0,'srate',20);
save_data(EEG_sim_test_case_1b,'test_case_1b')

%% Case 1c:  eventA: y~1 (with overlap, long data i.e. srate = 1000, datalength = 60*60)
intercept = struct();
intercept.eventname = 'eventA';
intercept.type = 'intercept';
intercept.predictorName = 'intercept';
intercept.overlap = 0.4;
intercept.effectsize = 2;

signals{1} = intercept;
%signals{1}.range = nan;

EEG_sim_test_case_1c = simulate_data(signals,'noise',0,'srate',1000,'datalength',60*60); % 60 min
save_data(EEG_sim_test_case_1c,'test_case_1c')

%% Case 1d:  eventA: y~1 (no overlap, unusual sampling rate i.e. srate = 23)
intercept = struct();
intercept.eventname = 'eventA';
intercept.type = 'intercept';
intercept.predictorName = 'intercept';
intercept.overlap = -1;
intercept.effectsize = 2;

signals{1} = intercept;
%signals{1}.range = nan;

EEG_sim_test_case_1d = simulate_data(signals,'noise',0,'srate',23);
save_data(EEG_sim_test_case_1d,'test_case_1d')

%% Case 2a: eventA: y~1+cat (without overlap)
intercept = struct();
intercept.eventname = 'eventA';
intercept.type = 'intercept';
intercept.predictorName = 'intercept';
intercept.overlap = -1;
intercept.effectsize = 2;

cat= struct();
cat.eventname = 'eventA';
cat.type = '1x2';
cat.predictorName = 'conditionA';
cat.overlap = -1;
cat.effectsize = 3;
%cat.range = [0,100];

signals{1} = intercept;
signals{1}(2) = cat;
%signals{1}.range = nan;

EEG_sim_test_case_2a = simulate_data(signals,'noise',0,'srate',20);
save_data(EEG_sim_test_case_2a,'test_case_2a')

%% Case 2b: eventA: y~1+cat (with biased overlap (i.e. depending on condition))
intercept = struct();
intercept.eventname = 'eventA';
intercept.type = 'intercept';
intercept.predictorName = 'intercept';
intercept.overlap = 0.4;
intercept.effectsize = 2;

cat= struct();
cat.eventname = 'eventA';
cat.type = '1x2';
cat.predictorName = 'conditionA';
cat.overlap = 0.2;
cat.effectsize = 3;
%cat.range = [0,100];

signals{1} = intercept;
signals{1}(2) = cat;
%signals{1}.range = nan;

EEG_sim_test_case_2b = simulate_data(signals,'noise',0,'srate',20);
save_data(EEG_sim_test_case_2b,'test_case_2b')

%% Case 3a: eventA: y~1+cat+cont (without overlap)
intercept = struct();
intercept.eventname = 'eventA';
intercept.type = 'intercept';
intercept.predictorName = 'intercept';
intercept.overlap = -1;
intercept.effectsize = 2;
intercept.range = nan;

cat= struct();
cat.eventname = 'eventA';
cat.type = '1x2';
cat.predictorName = 'conditionA';
cat.overlap = -1;
cat.effectsize = 3;
cat.range = nan;

cont = struct();
cont.eventname = 'eventA';
cont.type = 'continuous';
cont.predictorName = 'continuousA';
cont.overlap = -1;
cont.effectsize = 4;
cont.range = [-1,1];

signals{1} = intercept;
signals{1}(2) = cat;
signals{1}(3) = cont;

EEG_sim_test_case_3a = simulate_data(signals,'noise',0,'srate',20);
save_data(EEG_sim_test_case_3a,'test_case_3a')

%% Case 3b: eventA: y~1+cat+cont (with overlap (intercept only))
intercept = struct();
intercept.eventname = 'eventA';
intercept.type = 'intercept';
intercept.predictorName = 'intercept';
intercept.overlap = 0.4;
intercept.effectsize = 2;
intercept.range = nan;

cat= struct();
cat.eventname = 'eventA';
cat.type = '1x2';
cat.predictorName = 'conditionA';
cat.overlap = -1;
cat.effectsize = 3;
cat.range = nan;

cont = struct();
cont.eventname = 'eventA';
cont.type = 'continuous';
cont.predictorName = 'continuousA';
cont.overlap = -1;
cont.effectsize = 4;
cont.range = [-1,1];

signals{1} = intercept;
signals{1}(2) = cat;
signals{1}(3) = cont;

EEG_sim_test_case_3b = simulate_data(signals,'noise',0,'srate',20);
save_data(EEG_sim_test_case_3b,'test_case_3b')

%% Case 4a: eventA: y~1, eventB: y~1 (without overlap)
intercept1 = struct();
intercept1.eventname = 'eventA';
intercept1.type = 'intercept';
intercept1.predictorName = 'intercept';
intercept1.overlap = -1;
intercept1.effectsize = 2;

intercept2 = struct();
intercept2.eventname = 'eventB';
intercept2.type = 'intercept';
intercept2.predictorName = 'intercept';
intercept2.overlap = -1;
intercept2.effectsize = 3;

signals{1} = intercept1;
signals{2} = intercept2;
            
EEG_sim_test_case_4a = simulate_data(signals,'noise',0,'srate',20);
save_data(EEG_sim_test_case_4a,'test_case_4a')

%% Case 4b: eventA: y~1, eventB: y~1 (with overlap)
intercept1 = struct();
intercept1.eventname = 'eventA';
intercept1.type = 'intercept';
intercept1.predictorName = 'intercept';
intercept1.overlap = -1;
intercept1.effectsize = 2;

intercept2 = struct();
intercept2.eventname = 'eventB';
intercept2.type = 'intercept';
intercept2.predictorName = 'intercept';
intercept2.overlap = 0.4;
intercept2.effectsize = 3;

signals{1} = intercept1;
signals{2} = intercept2;
            
EEG_sim_test_case_4b = simulate_data(signals,'noise',0,'srate',20);
save_data(EEG_sim_test_case_4b,'test_case_4b')

%% Function to save data
function save_data(EEG_sim,file_name)
    % Check whether there is already a data folder, otherwise create one
    data_path = '../test/data_new_testcases/';
    if ~exist(data_path, 'dir')
       mkdir(data_path)
    end
    
    writematrix(EEG_sim.data,strcat(data_path,file_name,'_data.csv'));
    writetable(struct2table(EEG_sim.event),strcat(data_path,file_name,'_events.csv'));
end