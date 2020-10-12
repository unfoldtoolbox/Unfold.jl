ufresult_all = [];
ufresult_all_a = [];
badSubject = [];

for sub = 1:40
    %%
    filename = ['C:/Users/behinger/Downloads/P3_clean/' num2str(sub) '_P3_shifted_ds_reref_ucbip_hpfilt_ica_weighted_clean.set'];
    filename_csv= ['C:/Users/behinger/Downloads/P3_clean/' num2str(sub) '_P3_shifted_ds_reref_ucbip_hpfilt_ica_weighted_clean.set.csv'];
    try
        csv = readtable(filename_csv);
    catch e
        warning(e.message)
        badSubject = unique([badSubject sub]);
        continue
    end
    
    
    EEG= pop_loadset(filename);
    EEG = pop_eegfiltnew(EEG,0.5,[]);
    
    winrej = uf_continuousArtifactDetectASR(EEG);
    
    % rename to be used in EEG
    EEG.event = table2struct(csv);
    for e  = 1:length(EEG.event)
        EEG.event(e).type = EEG.event(e).eventtype;
        
        % manual interaction
        if EEG.event(e).trialtype == "target"
            EEG.event(e).rtDis = 0;
            EEG.event(e).rtTar = EEG.event(e).rt;
        elseif EEG.event(e).trialtype == "distractor"
            EEG.event(e).rtDis = EEG.event(e).rt;
            EEG.event(e).rtTar = 0;
        end
        
    end
    badIx = cellfun(@(x)isnan(x),{EEG.event.rt});
    EEG.event(badIx) = [];
    %EEG = pop_reref(EEG);
    %%
    for formula = {'y~1+cat(trialtype)'
            'y~1+cat(trialtype)+rt'
            'y~1+cat(trialtype)+spl(rt,4)'
            ...%'y~1+cat(trialtype)+spl(rtDis,4) + spl(rtTar,4)'
            }'
        %for k = 1:2
        EEG = uf_designmat(EEG,'eventtypes',{'stimulus','button'},'formula',...
            {formula{1},'y~1+cat(trialtype)'});
       
        EEG = uf_timeexpandDesignmat(EEG,'timelimits',[-1 1]);
        EEG = uf_continuousArtifactExclude(EEG,'winrej',winrej);
        
        EEG = uf_glmfit(EEG,'channel',21);
        
        EEGe = uf_epoch(EEG,'timelimits',[-1 1]);
        predictAt = {{'rt',[200 250 300 350 400 450]}};
        EEGe = uf_glmfit_nodc(EEGe,'method','matlab');
        
        ufresult = uf_condense(EEGe);
        ufresult_a = uf_addmarginal(uf_predictContinuous(ufresult,'predictAt',predictAt));
  
%         EEG = uf_glmfit(EEG,'channel',21,'method','glmnet','glmnetalpha',0);
%         EEGe = uf_glmfit_nodc(EEGe,'method','glmnet','glmnetalpha',0,'channel',21);
        
        filename = sprintf('sub-%i_formula-%s.mat',sub,formula{1});
        if ~exist(fullfile('local','p3'),'dir')
            mkdir(fullfile('local','p3'))
        end
        save(fullfile('local','p3',filename),'ufresult_a')
        
    end
    
end
% for bSet = betaSetName
%     ufresult_all.(bSet{1}) = ufresult_all.(bSet{1})(:,:,:,setdiff(1:40,badSubject));
%     ufresult_all_a.(bSet{1}) = ufresult_all_a.(bSet{1})(:,:,:,setdiff(1:40,badSubject));
% end
% %%
%%
tmp_fn_p3 = dir(fullfile('local','p3','*.mat'));
tmp_fn_p3 = {tmp_fn_p3.name};
fn_p3 = cellfun(@(x)strsplit(x,'_'),tmp_fn_p3,'UniformOutput',false);
fn_p3 = cell2table(cat(1,fn_p3{:}),'VariableNames',{'sub','formula'});

%fn_p3 = parse_column(fn_p3,'overlap');
%fn_p3 = parse_column(fn_p3,'noise');
fn_p3.filename = tmp_fn_p3';
fn_p3.folder = repmat({'p3'},1,height(fn_p3))';
%%
all_b = nan(height(fn_p3),31,512,10);
all_bnodc = nan(height(fn_p3),31,512,10);
for r = 1:height(fn_p3)
    fprintf("Loading :%i/%i\n",r,height(fn_p3))
    
    tmp = load(fullfile('local',fn_p3.folder{r},fn_p3.filename{r}));
    b = tmp.ufresult_a.beta(:,:,:);
    b_nodc = tmp.ufresult_a.beta_nodc(:,:,:);
    if strcmp(fn_p3{r,'formula'},'formula-y~1+cat(trialtype).mat')
        b(:,:,8:10) = b(:,:,2:4);
        b(:,:,2:7) = repmat(b(:,:,1),1,1,6);
        b_nodc(:,:,8:10) = b_nodc(:,:,2:4);
        b_nodc(:,:,2:7) = repmat(b_nodc(:,:,1),1,1,6);
    end
    
    
    all_b(r,:,:,:) = b;
    all_bnodc(r,:,:,:) =  b_nodc;

    %fn_p3{r,'ufresult'} = tmp.ufresult_marginal;
end
fn_p3.beta = squeeze(all_b);
fn_p3.beta_nodc = squeeze(all_bnodc);

%%
plot_result(fn_p3(1,:),'channel',21)

%% GA
groupIx = findgroups(fn_p3.formula);
GA = splitapply(@(x)trimmean(x,0.2),fn_p3.beta,groupIx);
GA_nodc = splitapply(@(x)trimmean(x,0.2),fn_p3.beta_nodc,groupIx);
fn_p3_ga = table(unique(fn_p3.formula),GA,GA_nodc,'VariableNames',{'formula','beta','beta_nodc'});
fn_p3_ga.folder = repmat({'p3'},1,height(fn_p3_ga))';
fn_p3_ga.filename = fn_p3{1:3,'filename'}
plot_result(fn_p3_ga(2,:),'channel',21)


%%
%uf_plotParam(ufresult_a,'channel',21)