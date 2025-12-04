T = readtable("parametric_total.csv");

geom = T(:, {'a_mm','b_mm','c_mm','t_mm','n1','n2'});

[uniqueGeom, ia2, ic2] = unique(geom, 'rows');

duplicate_geom_idx = setdiff(1:height(T), ia2);

if isempty(duplicate_geom_idx)
    disp("No duplicate geometries found.");
else
    disp("Duplicate geometries found at these rows:");
    disp(duplicate_geom_idx');

    disp("Duplicate geometry parameter rows:");
    disp(geom(duplicate_geom_idx, :));
end
T_clean = T(ia2, :);   % keep only unique geometries
writetable(T_clean, "cleaned_file.csv");

disp("Cleaned file saved as cleaned_file.csv");