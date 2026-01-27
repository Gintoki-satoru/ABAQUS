data = readmatrix('S12_data.xlsx');

% Sphere - [0,-30,60,30]
x = data(2:end, 1);
y = data(2:end, 2);

% Plot
figure
plot(x, y, 'o')
grid on
xlim([0 250]);
ylim([-5 1]);
xlabel('X')
ylabel('S12')
title('S12 vs X - Sphere - [0,-30,60,30]')

%% Sphere - [0,45,-45,0,0,-45,45,0]
data = readmatrix('S12_data.xlsx');
x = data(2:end, 3);
y = data(2:end, 4);

% Plot
figure
plot(x, y, 'o')
grid on
xlim([40 250]);
xlabel('X')
ylabel('S12')
title('S12 vs X - Sphere - [0,45,-45,0,0,-45,45,0]')

%% Sphere - [90,30,-30,90]
data = readmatrix('S12_data.xlsx');
x = data(2:end, 9);
y = data(2:end, 10);

% Plot
figure
plot(x, y, 'o')
grid on
xlim([0 250]);
xlabel('X')
ylabel('S12')
title('S12 vs X - Sphere - [90,30,-30,90]')

%% Cigar - [0,-30,60,30]
data = readmatrix('S12_data.xlsx');
x = data(2:end, 5);
y = data(2:end, 6);

% Plot
figure
plot(x, y, 'o')
grid on
xlim([0 116]);
xlabel('X')
ylabel('S12')
title('S12 vs X - Cigar - [0,-30,60,30]')

%% Cigar - [0,45,-45,0,0,-45,45,0]
data = readmatrix('S12_data.xlsx');
x = data(2:end, 7);
y = data(2:end, 8);

% Plot
figure
plot(x, y, 'o')
grid on
xlim([0 116]);
xlabel('X')
ylabel('S12')
title('S12 vs X - Cigar - [0,45,-45,0,0,-45,45,0]')