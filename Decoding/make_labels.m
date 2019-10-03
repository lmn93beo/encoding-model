function [binned_labels, valid_items] = make_labels(data, decoding_type, ncells, validTrials)
% Make labels for decoding
ntrials = numel(data.response.choice);


if strcmp(decoding_type, 'stim')
    % L/R stim
    ltrials = find(data.stimuli.loc(1:ntrials) == 2);
    rtrials = find(data.stimuli.loc(1:ntrials) == 1);

    lvalid = intersect(ltrials, validTrials);
    rvalid = intersect(rtrials, validTrials);

    valid_items = union(lvalid, rvalid);
    binned_labels.side = cell(1, ncells);

    for i = 1:ncells
        binned_labels.side{i} = data.stimuli.loc(valid_items);
    end
    
elseif strcmp(decoding_type, 'choice')
    % L/R stim
    ltrials = find(data.response.choice(1:ntrials) == 2);
    rtrials = find(data.response.choice(1:ntrials) == 1);

    lvalid = intersect(ltrials, validTrials);
    rvalid = intersect(rtrials, validTrials);

    valid_items = union(lvalid, rvalid);
    binned_labels.side = cell(1, ncells);

    for i = 1:ncells
        binned_labels.side{i} = data.response.choice(valid_items);
    end
    
elseif strcmp(decoding_type, 'reward')
    rewards = find(data.response.reward > 0);
    norew = find(data.response.reward == 0 & data.response.choice ~= 5);
    reward_valid = intersect(rewards, validTrials);
    norew_valid = intersect(norew, validTrials);
    
    valid_items = union(reward_valid, norew_valid);
    
    binned_labels.side = cell(1, ncells);

    for i = 1:ncells
        binned_labels.side{i} = data.response.reward(valid_items);
    end
end