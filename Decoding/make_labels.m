function [binned_labels, valid_items] = make_labels(data, decoding_type, ncells, validTrials)
% Make labels for decoding
ntrials = numel(data.response.choice);

switch decoding_type
    case 'stim'    % L/R stim
        ltrials = find(data.stimuli.loc(1:ntrials) == 2);
        rtrials = find(data.stimuli.loc(1:ntrials) == 1);
        
        [binned_labels, valid_items] = make_labels_helper(data.stimuli.loc, ltrials, rtrials, ncells, validTrials);
    
    case 'choice'
        ltrials = find(data.response.choice(1:ntrials) == 2);
        rtrials = find(data.response.choice(1:ntrials) == 1);

        [binned_labels, valid_items] = make_labels_helper(data.response.choice, ltrials, rtrials, ncells, validTrials);
    
    case 'reward'
        rewards = find(data.response.reward > 0);
        norew = find(data.response.reward == 0 & data.response.choice ~= 5);
        
        [binned_labels, valid_items] = make_labels_helper(data.response.reward, rewards, norew, ncells, validTrials);
        
    case 'left_difficulty'
        locs = data.stimuli.loc(1:ntrials);
        opp = data.stimuli.opp_contrast(1:ntrials);
        left_corr_easy = find(data.response.choice == 2 & locs == 2 & opp < 0.2);
        left_corr_difficult = find(data.response.choice == 2 & locs == 2 & opp > 0.2);
        [binned_labels, valid_items] = make_labels_helper(opp > 0.2, left_corr_easy, left_corr_difficult, ncells, validTrials);
        
    case 'right_difficulty'
        locs = data.stimuli.loc(1:ntrials);
        opp = data.stimuli.opp_contrast(1:ntrials);
        left_corr_easy = find(data.response.choice == 1 & locs == 1 & opp < 0.2);
        left_corr_difficult = find(data.response.choice == 1 & locs == 1 & opp > 0.2);
        [binned_labels, valid_items] = make_labels_helper(opp > 0.2, left_corr_easy, left_corr_difficult, ncells, validTrials);
end
end

function [binned_labels, valid_items] = make_labels_helper(series, Atrials, Btrials, ncells, validTrials)
    Avalid = intersect(Atrials, validTrials);
    Bvalid = intersect(Btrials, validTrials);

    valid_items = union(Avalid, Bvalid);
    binned_labels.side = cell(1, ncells);

    for i = 1:ncells
        binned_labels.side{i} = series(valid_items);
    end

end
