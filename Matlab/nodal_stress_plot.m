% Ply 2&3
data = readmatrix('Node stress.xlsx');
x = data(2:end, 1);
y = data(2:end, 3);
z = data(2:end, 4); 

% Plot
figure
plot(x, y, 'o')
grid on
% xlim([0 250]);
% ylim([-5 1]);
xlabel('X')
ylabel('S12')
title('S12 vs X - Sphere - [90,30,-30,90] - Ply 2&3')

figure
plot(x, z, 'o')
grid on
% xlim([0 250]);
% ylim([-5 1]);
xlabel('X')
ylabel('S22')
title('S22 vs X - Sphere - [90,30,-30,90] - Ply 2&3')

%% Ply 3&4
data = readmatrix('Node stress.xlsx');

x = data(2:end, 5);
y = data(2:end, 7);
z = data(2:end, 8); 

% Plot
figure
plot(x, y, 'o')
grid on
% xlim([0 250]);
% ylim([-5 1]);
xlabel('X')
ylabel('S12')
title('S12 vs X - Sphere - [90,30,-30,90] - Ply 3&4')

figure
plot(x, z, 'o')
grid on
% xlim([0 250]);
% ylim([-5 1]);
xlabel('X')
ylabel('S22')
title('S22 vs X - Sphere - [90,30,-30,90] - Ply 3&4')

%% Ply 1&2
data = readmatrix('Node stress.xlsx');

x = data(2:end, 9);
y = data(2:end, 11);
z = data(2:end, 12); 

% Plot
figure
plot(x, y, 'o')
grid on
% xlim([0 250]);
% ylim([-5 1]);
xlabel('X')
ylabel('S12')
title('S12 vs X - Sphere - [90,30,-30,90] - Ply 1&2')

figure
plot(x, z, 'o')
grid on
% xlim([0 250]);
% ylim([-5 1]);
xlabel('X')
ylabel('S22')
title('S22 vs X - Sphere - [90,30,-30,90] - Ply 1&2')

