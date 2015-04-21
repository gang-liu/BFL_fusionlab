SBJ_CELL = {'IBSR_01' 'IBSR_03' 'IBSR_05' ...
        'IBSR_07' 'IBSR_09' 'IBSR_11' 'IBSR_13' ...
        'IBSR_15' 'IBSR_17' 'IBSR_02' ...
        'IBSR_04' 'IBSR_06' 'IBSR_08' ...
        'IBSR_10' 'IBSR_12' 'IBSR_14' ...
        'IBSR_16' 'IBSR_18'...
        };
addpath('../misc');
DATA_DIR = '../data/';
Output_dir = '../temp/';

labels = [17 53];

dice_stack = nan(10, 4);
for sub = 1:10
    sub
    sbj = SBJ_CELL{sub};
    man_seg = MRIread([DATA_DIR '/' sbj '_man_seg.mgz']);
    tmp_cnt = 1;
    for sig = [1 5 10 50]
        
        aseg = MRIread([Output_dir '/' sbj '_mri_norm.mgz_BFL_seg_sig_lf_' num2str(sig) '.mgz']);
        
        tmp_dice = computeLabelDice(man_seg.vol, aseg.vol, labels)
        dice_stack(sub, tmp_cnt) = mean(tmp_dice);
        tmp_cnt = tmp_cnt + 1;
    end
        
    
end