for sub = 1:40
%%
    filename = ['C:/Users/behinger/Downloads/P3_clean/' num2str(sub) '_P3_shifted_ds_reref_ucbip_hpfilt_ica_weighted_clean.set'];
    filename_csv= ['C:/Users/behinger/Downloads/P3_clean/' num2str(sub) '_P3_shifted_ds_reref_ucbip_hpfilt_ica_weighted_clean.set.csv'];
    csv = readtable(filename_csv);

    EEG= pop_loadset(filename);
    EEG.event = table2struct(csv);
    for e  = 1:length(EEG.event)
    EEG.event(e).type = EEG.event(e).eventtype;
    end
    %EEG = pop_reref(EEG);
    EEG = uf_designmat(EEG,'eventtypes',{'stimulus','button'},'formula',{'y~1+cat(trialtype)','y~1+cat(trialtype)'});
    EEG = uf_timeexpandDesignmat(EEG,'timelimits',[-1 1]);
    EEG = uf_glmfit(EEG,'method','pinv');
    EEGe = uf_epoch(EEG,'timelimits',[-1 1]);
    EEGe = uf_glmfit_nodc(EEGe);
    ufresult = uf_condense(EEGe);
end