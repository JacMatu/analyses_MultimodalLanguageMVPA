function ds = setTargetsChunksLabels(opt, ds, stim, modalityNb, nbRun)

    % set chunks (runs by trial_type), target (stimulation type - modality),
    % names (stimulation type name)
    trialsPerRun = length(stim) * opt.cosmomvpa.nbTrialRepetition;

    chunks = repmat((1:nbRun)', 1, opt.cosmomvpa.nbTrialRepetition)';
    chunks = chunks(:);
    chunks = repmat(chunks, trialsPerRun, 1);

    targets = repmat(stim', 1, nbRun)';
    targets = targets(:);


    modalityNb =  repmat(modalityNb,(nbRun*opt.cosmomvpa.nbTrialRepetition),1);
    modalityNb = modalityNb(:);

    ds.sa.targets = targets;
    ds.sa.chunks = chunks;
    ds.sa.modality = modalityNb;
    
%     horzcat(ds.sa.chunks, ds.sa.targets, ds.sa.modality)

