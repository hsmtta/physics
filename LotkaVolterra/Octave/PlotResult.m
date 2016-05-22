%
% (c) 2016 Takahiro Hashimoto
%

data = dlmread('../result.dat', '', 2, 0);

figure(1);
plot(data(:,1), data(:,2), data(:,1), data(:,3));
grid on;
xlabel('Time');
ylabel('Population');

figure(2);
plot(data(:,1), data(:, 4));
grid on;
xlabel('Time');
ylabel('V');
axis([0 100 -3 0]);
