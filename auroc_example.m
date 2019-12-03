xtest = rand(1, 100);
ytest = rand(1, 100) + 3;

mix = [xtest ytest];
labels = [ones(1, 100) ones(1, 100) * 2];

[a, b, c, auroc] = perfcurve(labels, mix, 2);