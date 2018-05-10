fprintf('level=%d, ', level);
fprintf('time spent on K: %6.1e seconds, ',time_stiffness_matrix);
if exist('time_mass_matrix','var')
    fprintf('time spent on M: %6.1e seconds, ',time_mass_matrix);
end
fprintf('size of square matrix = %d ',row);
fprintf('\n');