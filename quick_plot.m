dates = {'20190724', '20190725', '20190726', '20190727',...
    '20190730', '20190801', '20190802', '20190806', '20190807'};

depths = [367.3, 338.3, 285.5, 373.15, 276.4, 322.7, 345.35, 382.75, 275.83];
n1 = [254, 110, 98, 261, 98, 76, 210, 155, 127];
n2 = [24, 32, 47, 35, 29, 3, 67, 28, 19];
n3 = [20, 3, 17, 13, 6, 23, 14, 10, 6];
total = n1 + n2 + n3;
plot(n1 ./ total, 'o');
hold on;
plot(n2 ./ total, 'o');
plot(n3 ./ total, 'o');


% acc = [77.29, 91.5, 86.5, 68.4, 83.5, 82.4, 79.4];
% t = [0.983, 3.197, 2.377, 1.229, 3.606, 0.573, 1.39];
% 
% acc = [90, 70, 77, 90, 60, 63, 73, 53]; 
% 
% plot(acc, 'o')