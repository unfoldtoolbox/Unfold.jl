for k = 2;[1 2 3 15]
    
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
            
            if 1 == 0
                %%
                EEG = uf_designmat(EEGsim,'formula','y~1+conditionA+conditionB','eventtypes','stimulus2');
                tic
                EEG = uf_timeexpandDesignmat(EEG,'timelimits',[-1,1]);
                toc
                tic
                EEG = uf_glmfit(EEG);
                toc
            end
        case 15
            
            EEGsim = simulate_test_case(15,'noise',0,'basis','posneg');
    end
    dlmwrite(sprintf('../test/data/testCase%i_data.csv',k),EEGsim.data,'precision','%.10f')
    writetable(struct2table(EEGsim.event),sprintf('../test/data/testCase%i_events.csv',k))
    
    
    
end