%% --- SETTINGS ---
numSamples = 50000;
numBest    = 15;

% Parameter ranges
range_a  = [100 1000];
range_b  = [100 1000];
range_c  = [100 1000];
range_n1 = [0.1 1];
range_n2 = [0.1 1];

%% --- STORAGE ---
results = [];

%% --- RANDOM SEARCH OPTIMIZATION ---
for i = 1:numSamples
    % Random parameters
    a  = rand_range(range_a);
    b  = rand_range(range_b);
    c  = rand_range(range_c);
    n1 = rand_range(range_n1);
    n2 = rand_range(range_n2);

    try
        V = superellipsoid_volume(a, b, c, n1, n2);
    catch
        % Skip invalid gamma combinations
        continue
    end

    % η = volume / bounding-box volume = V / (2a*2b*2c)
    eta = V / (8 * a * b * c);

    % Physically valid only for 0 < η ≤ 1
    if eta <= 0 || eta > 1
        continue
    end

    % Store result
    results = [results; eta, a, b, c, n1, n2];
end

%% --- SORT AND SELECT BEST ---
% Sort by efficiency descending
results = sortrows(results, 1, "descend");

% Keep only highest-performing unique structures
uniqueResults = unique(round(results, 4), 'rows', 'stable'); 
finalResults  = uniqueResults(1:min(numBest, size(uniqueResults,1)), :);

%% --- DISPLAY RESULTS ---
disp("Top high-packing-efficiency super-ellipsoid shapes:");
disp("   η        a        b        c        n1        n2");
disp(finalResults);


%% --- HELPER FUNCTION ---
function x = rand_range(range)
    x = range(1) + rand() * (range(2) - range(1));
end

function V = superellipsoid_volume(a, b, c, n1, n2)
    e1 =2/n1;
    e2 =2/n2;
    V = 8 * a * b * c * (gamma(1 + 1/e1)^2 * gamma(1 + 1/e2)) / ...
        (gamma(1 + 2/e1) * gamma(1 + (1/e2 + 2/e1)));
end
