% =========================================================================
%
%  This program triggers an assembly test for an elastic body. The
%  assembly time of the stiffness matrix is stored depending on a choice of
%  the finite elements and mesh density. 
%
% ======================================================================
%

% number of mesh levels 
levels=1;

% P1 elements
fprintf('P1 elements: \n')
level_time_P1=zeros(1,levels);
level_size_P1=zeros(1,levels);
for level=0:levels
    [level_time_P1(level+1),level_size_P1(level+1)]=elasticity_fem('P1',level,level==0);
end
fprintf('\n')

% P2 elements
fprintf('P2 elements: \n')
level_time_P2=zeros(1,levels);
level_size_P2=zeros(1,levels);
for level=0:levels
    [level_time_P2(level+1),level_size_P2(level+1)]=elasticity_fem('P2',level,level==0);
end
fprintf('\n')

% Q1 elements
fprintf('Q1 elements: \n')
level_time_Q1=zeros(1,levels);
level_size_Q1=zeros(1,levels);
for level=0:levels
    [level_time_Q1(level+1),level_size_Q1(level+1)]=elasticity_fem('Q1',level,level==0);
end
fprintf('\n')

% Q2 elements
fprintf('Q2 elements: \n')
level_time_Q2=zeros(1,levels);
level_size_Q2=zeros(1,levels);
for level=0:levels
    [level_time_Q2(level+1),level_size_Q2(level+1)]=elasticity_fem('Q2',level,level==0);
end

% output information
fprintf('\n')
for level=0:levels
    fprintf('%d ', level);
    fprintf('& ');
    fprintf('%d ', level_size_P1(level+1));
    fprintf('& ');
    fprintf('%d ', level_size_P2(level+1));
    fprintf('& ');
    fprintf('%d ', level_size_Q1(level+1));
    fprintf('& ');
    fprintf('%d ', level_size_Q2(level+1));
    fprintf('& ');
    fprintf('%2.2f ', level_time_P1(level+1));
    fprintf('& ');
    fprintf('%2.2f ', level_time_P2(level+1));
    fprintf('& ');
    fprintf('%2.2f ', level_time_Q1(level+1));
    fprintf('& ');
    fprintf('%2.2f ', level_time_Q2(level+1));
    fprintf('\\\\');
    fprintf('\n');
end