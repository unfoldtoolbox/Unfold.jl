% Fuction to simulate test cases and export them
% you need to have the matlab-Unfold toolbox in your path & initialized
% Testcase 3 is multiSubject
% Testcase 4 is long without overlap (srate = 1000)
% Testcase 5 is long with overlap (srate = 1000)
% Testacse 1 is without overlap
% testcase 2 should have overlap
% Tastcase 15 is complex
for k = 6;5;[1 2 3 4 5 15];

    rng(1);

    switch k

        case 1
            %%

            % Basis functions
            intercept = struct();
            intercept.eventname = 'stimulus2';
            intercept.type = 'intercept';
            intercept.predictorName = 'intercept';
            intercept.overlap = 0; % note that 0, means that the AVERAGE overlap is exactly one signal-length (i.e. 1s in our simulations). That means if we do not want ANY overlap, we specify -1 here!
            intercept.effectsize = 3;

            cont = struct();
            cont.eventname = 'stimulusA';
            cont.type = 'continuous';
            cont.overlap = 0;
            cont.predictorName = 'continuousA';
            cont.effectsize = 1;
            cont.range = [0,100];

            signals{1} = intercept;
            signals{1}.range = nan;
            signals{1}.overlap = -1;
            signals{1}(2) = cont;
            signals{1}(2).overlap = -1;
            signals{1}(2).effectsize= 2.5;
            signals{1}(3) = cont;
            signals{1}(3).type = "continuous";
            signals{1}(3).predictorName = 'continuousB';
            signals{1}(3).overlap = -1;
            signals{1}(3).effectsize= -1.5;

            EEGsim = simulate_data(signals,'noise',0);

            %' Formula: y ~ 1  + continuousA + continuousB'
            %' Overlap: y ~ -1 +  -1         + -1 => no overlap
            %' EffectS: y ~ 3  +  2.5        + -1.5
            %' Shape: 'box', noise: 0

            %' Overlap: y ~ 0.5 +  0        + 0 => overlap
            %'
            % stim: y~1
            % button: y~1

            % event: y~1 + cat(is_button)
            % event (stim) => y~stimERP
            % event (button) => y~stimERP + (stimERP-buttonEffect)
        case 2
            %%
            % Basis functions
            intercept = struct();
            intercept.eventname = 'stimulus2';
            intercept.type = 'intercept';
            intercept.predictorName = 'intercept';
            intercept.overlap = 0.5; % note that 0, means that the AVERAGE overlap is exactly one signal-length (i.e. 1s in our simulations). That means if we do not want ANY overlap, we specify -1 here!
            intercept.effectsize = 3;

            cat= struct();
            cat.eventname = 'stimulusA';
            cat.type = '1x2';
            cat.overlap = 0;
            cat.predictorName = 'conditionA';
            cat.effectsize = 1;
            cat.range = [0,100];

            signals{1} = intercept;
            signals{1}.range = nan;
            signals{1}.overlap = 0;
            signals{1}(2) = cat;
            signals{1}(2).overlap = 0;
            signals{1}(2).effectsize= 2.5;
            signals{1}(3) = cat;
            signals{1}(3).predictorName = 'conditionB';
            signals{1}(3).overlap = 0.1;
            signals{1}(3).effectsize= -1.5;

            EEGsim = simulate_data(signals,'noise',0.5,'basis','posneg');

            %' Formula: y ~ 1  + conditionA + conditionB'
            %' Overlap: y ~ 0  +  0         + 0.1
            %' EffectS: y ~ 3  +  2.5       + -1.5
            %' Shape: 'posneg', noise: 0.5
        case 3
            %%
            %%
            rng(1)
            input = [];
            for subj = 1:50
                % small noise
                input{subj} = simulate_data_lmm_v2('noise',0,'u_noise',0,'noise_components',0,...
                    'srate',10,'n_events',20+2*subj,...
                    'b_p1_2x2',[0,0,0,0],... % P1: Intercept, MainA, MainB, Inter - effect coded beta
                    'u_p1_2x2',[1,0,0,0],... %   P1: Subject variability
                    ...'b_p1_2x2',[10,5,10,0],... % P1: Intercept, MainA, MainB, Inter - effect coded beta
                    ...'u_p1_2x2',[3 ,3,3,0],... %   P1: Subject variability
                    ...
                    'b_n1_2x2',[0,0,0,0],...
                    'u_n1_2x2',[0 ,0 ,0,0],...
                    ...
                    'b_p3_2x2',[0,0,0,0],...
                    'u_p3_2x2',[0 ,0,0,0],...
                    ...
                    'u_p1_item',[0], ... % Item effect strength
                    'u_n1_item',[0], ...
                    'u_p3_item',[0], ...
                    'simulationtype','ideal', ...
                    'overlaptype','uniform',...
                    'overlapparam',[1 0 0 0;2 0 0 0],...
                    'randomItem',0 ...
                    );


            end



            %%

            cfgDesign = struct();
            cfgDesign.inputGroupingName='subject';
            cfgDesign.eventtypes= 'sim';
            %             cfgDesign.formula = 'y~1+cat(condA)+cat(condB)+(1+condA+condB|subject)';%+(1|stimulus)';
            cfgDesign.formula = 'y~1+cat(condA)+cat(condB)+(1|subject)';%+(1|stimulus)';
            cfgDesign.codingschema = 'effects';

            [EEG,EEG_fixef] = um_designmat(input,cfgDesign);

            EEG= um_timeexpandDesignmat(EEG,'timelimits',[-.15,0.4]);
            EEG = um_mmfit(EEG,input)
            x = um_condense(EEG)
            %%
            EEGsim.data = cellfun(@(x)squeeze(x.data(1,:,:)),input,'UniformOutput',0);

            y = cat(2,EEGsim.data{:});
            EEGsim.event = EEG_fixef.event;


        case 4
            % long example
            %%
            % Basis functions
            intercept = struct();
            intercept.eventname = 'stimulus2';
            intercept.type = 'intercept';
            intercept.predictorName = 'intercept';
            intercept.overlap = 0.5; % note that 0, means that the AVERAGE overlap is exactly one signal-length (i.e. 1s in our simulations). That means if we do not want ANY overlap, we specify -1 here!
            intercept.effectsize = 3;

            cat= struct();
            cat.eventname = 'stimulusA';
            cat.type = '1x2';
            cat.overlap = 0.5;
            cat.predictorName = 'conditionA';
            cat.effectsize = 1;
            cat.range = [0,100];

            signals{1} = intercept;
            signals{1}.range = nan;
            signals{1}.overlap = -1;
            signals{1}(2) = cat;
            signals{1}(2).overlap = 0.5;
            signals{1}(2).effectsize= 2.5;
            signals{1}(3) = cat;
            signals{1}(3).predictorName = 'conditionB';
            signals{1}(3).overlap = 0.1;
            signals{1}(3).effectsize= -1.5;

            EEGsim = simulate_data(signals,'noise',0.5,'basis','posneg','srate',1000,'datalength',50*60);% 50min

            %' Formula: y ~ 1  + conditionA + conditionB'
            %' Overlap: y ~ -1 +  0.5       + 0.1
            %' EffectS: y ~ 3  +  2.5       + -1.5
            %' Shape: 'posneg', noise: 0.5

            if 1 == 0
                %%
                EEG = uf_designmat(EEGsim,'formula','y~1+conditionA+conditionB','eventtypes','stimulus2');
                tic
                EEG = uf_timeexpandDesignmat(EEG,'timelimits',[-1,1]);
                toc
                tic
                EEG = uf_glmfit(EEG);
                EEG = uf_epoch(EEG,'timelimits',[-1,1])
                EEG = uf_glmfit_nodc(EEG)
                uf_plotParam(uf_condense(EEG))
                toc
            end
        case 5
             % long example
            %%
            % Basis functions
            intercept = struct();
            intercept.eventname = 'stimulus2';
            intercept.type = 'intercept';
            intercept.predictorName = 'intercept';
            intercept.overlap = 0.5; % note that 0, means that the AVERAGE overlap is exactly one signal-length (i.e. 1s in our simulations). That means if we do not want ANY overlap, we specify -1 here!
            intercept.effectsize = 3;

            cat= struct();
            cat.eventname = 'stimulusA';
            cat.type = '1x2';
            cat.overlap = 0.5;
            cat.predictorName = 'conditionA';
            cat.effectsize = 1;
            cat.range = [0,100];

            signals{1} = intercept;
            signals{1}.range = nan;
            signals{1}.overlap = 0.5;
            signals{1}(2) = cat;
            signals{1}(2).overlap = 0.5;
            signals{1}(2).effectsize= 2.5;
            signals{1}(3) = cat;
            signals{1}(3).predictorName = 'conditionB';
            signals{1}(3).overlap = 0.1;
            signals{1}(3).effectsize= -1.5;

            EEGsim = simulate_data(signals,'noise',0.5,'basis','posneg','srate',1000,'datalength',50*60);% 50min

            %' Formula: y ~ 1  + conditionA + conditionB'
            %' Overlap: y ~ 0.5 +  0.5       + 0.1
            %' EffectS: y ~ 3  +  2.5       + -1.5
            %' Shape: 'posneg', noise: 0.5

            if 1 == 0
                %%
                EEG = uf_designmat(EEGsim,'formula','y~1+conditionA+conditionB','eventtypes','stimulus2');
                tic
                EEG = uf_timeexpandDesignmat(EEG,'timelimits',[-1,1]);
                toc
                tic
                EEG = uf_glmfit(EEG);
                EEG = uf_epoch(EEG,'timelimits',[-1,1])
                EEG = uf_glmfit_nodc(EEG)
                uf_plotParam(uf_condense(EEG))
                toc
            end
        case 15

            EEGsim = simulate_test_case(15,'noise',0,'basis','posneg');
            % Stim1 : y~1
            % stim2 : y~1+cat+cont
            % stim3 : y~1+splA+splB

        case 6 %'lmm_sub-item_realistic'
                        rng(1)
            input = [];
            for subj = 1:20
                % small noise
                input{subj} = simulate_data_lmm_v2('noise',0,'u_noise',0,'noise_components',0,...
                    'srate',100,'n_events',20+2*subj,...
                    'b_p1_2x2',[1,0.3,0,0],... % P1: Intercept, MainA, MainB, Inter - effect coded beta
                    'u_p1_2x2',[1,0.1,0,0],... %   P1: Subject variability
                    ...'b_p1_2x2',[10,5,10,0],... % P1: Intercept, MainA, MainB, Inter - effect coded beta
                    ...'u_p1_2x2',[3 ,3,3,0],... %   P1: Subject variability
                    ...
                    'b_n1_2x2',[-0.2,1,0,0],...
                    'u_n1_2x2',[1 ,1 ,0,0],...
                    ...
                    'b_p3_2x2',[3,0,0,0],...
                    'u_p3_2x2',[2,2,0,0],...
                    ...
                    'u_p1_item',[2], ... % Item effect strength
                    'u_n1_item',[3], ...
                    'u_p3_item',[0], ...
                    'n_items',[20], ...
                    'simulationtype','ideal_hanning', ...
                    'overlaptype','uniform',...
                    'overlapparam',[1 0 0 0;2 0 0 0],...
                    'overlapminimum',1,... # deactivating overlap!
                    'randomItem',1 ...
                    );
            end
            cfgDesign = struct();
            cfgDesign.inputGroupingName='subject';
            cfgDesign.eventtypes= 'sim';
            %             cfgDesign.formula = 'y~1+cat(condA)+cat(condB)+(1+condA+condB|subject)';%+(1|stimulus)';
            cfgDesign.formula = 'y~1+cat(condA)+(1+condA|subject)+(1+condA|stimulus)';
            cfgDesign.codingschema = 'effects';

            [EEG,EEG_fixef] = um_designmat(input,cfgDesign);

%             EEG= um_timeexpandDesignmat(EEG,'timelimits',[-.15,0.4]);
%             EEG = um_mmfit(EEG,input)
%             x = um_condense(EEG)
            %
            EEGsim.data = cellfun(@(x)squeeze(x.data(1,:,:)),input,'UniformOutput',0);

            y = cat(2,EEGsim.data{:});
            EEGsim.event = EEG_fixef.event;

    end
    dlmwrite(sprintf('../test/data/testCase%i_data.csv',k),EEGsim.data,'precision','%.10f')
    writetable(struct2table(EEGsim.event),sprintf('../test/data/testCase%i_events.csv',k))



end
