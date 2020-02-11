for k = [1 15]
    
    
    switch k
        case 1
            
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
            signals{1}(3).type = "1x2";
            signals{1}(3).predictorName = 'conditionA';
            signals{1}(3).overlap = -1;
            signals{1}(3).effectsize= -1.5;
            
            EEGsim = simulate_data(signals,'noise',0);
            
        case 15
            
            EEGsim = simulate_test_case(15,'noise',0,'basis','box');
    end
    csvwrite(sprintf('../test/testCase%i_data.csv',k),EEGsim.data)
    writetable(struct2table(EEGsim.event),sprintf('../test/testCase%i_events.csv',k))
    
    
    
end