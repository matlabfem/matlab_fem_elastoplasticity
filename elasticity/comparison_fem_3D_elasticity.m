clear all;
close all;

% add_paths; 
addpath(genpath('3rd_party_elasticity_codes_for_testing'),'elasticity_3D');

fprintf('Benchmark from the article of J. Alberty, C. Carstensen, S. Funken, R. Klose: \n')
figure(1); fem_lame3d;   %running the modified code of Alberty, Carstensen, Funken, Klose

fprintf('\n')

fprintf('Benchmark from the article of J. Koko: \n')
figure(2); demo_elas;    %running the modified code of Koko

rmpath(genpath('3rd_party_elasticity_codes_for_testing'),'elasticity_3D');



