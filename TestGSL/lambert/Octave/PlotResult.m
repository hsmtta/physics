%
% (c) 2016 Takahiro Hashimoto
%

data = dlmread('../build/result.dat', '', 1, 0);

plot(data(:,1), data(:,2));
grid on;
xlabel('x');
ylabel('W(x)');
