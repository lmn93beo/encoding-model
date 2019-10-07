load('TB146_decoding\TB146_choice_raw_20190727.mat');

labels = binned_labels.side{1};
left = find(labels == 2);
right = find(labels == 1);

% Plot left and right separately
cellid = 1;

cellarr = binned_data{cellid};
figure;
subplot('121');
plot(cellarr(left, :)');
subplot('122');
plot(cellarr(right, :)');
