%
% Linear elasticity  2D/3D
% 
% 2D: \Omega=(0 35)*(-5 5)
%     u=0 on x=0 
%     f=(0,-.75) on Omega
%     g=(0,10) on x=35
% 3D: \Omega=(-5 5)*(0 35)*(-5 5)
%     u=0 on y=0
%     f=(0,0,-.75) volume forces
%     g=(0,0,10) on y=35
%----------------------------------------------------------
E=30000; nu=0.4;
d=3;
fz=-0.75;
gz=10;
 
% mesh generation
if (d==2)
    ax=0; bx=35; ay=-5; by=5;  
    nx=36; ny=11;
    [p,t,bpx,bpy]=kpde2dumsh(ax,bx,ay,by,nx,ny);
    % Dirichlet at x=0
    ibcd=[2*bpy(:,1)-1; 2*bpy(:,1)];
    % Neumann at x=35
    ibn=bpy(:,2); y=p(ibn,2); [~,iy]=sort(y); ibn=ibn(iy); ne=length(ibn);
    eneum=[ibn(1:ne-1) ibn(2:ne)];
    ibcn=ibn;
    clear bpx bpy ibn y iy
elseif (d==3)
    ax=-5; bx=5; ay=0; by=35; az=-5; bz=5;
    nx=11; ny=26; nz=11;
    %nx=nx*2; ny=ny*2; nz=nz*2;                                              %added by Jan Valdman to make the mesh finer! 
    [p,t,bpx,bpy,bpz]=kpde3dumsh(ax,bx,ay,by,az,bz,nx,ny,nz);
    % Dirichlet at y=0
    ibcd=[3*bpy(:,1)-2; 3*bpy(:,1)-1; 3*bpy(:,1)];
    % Neumann condition at y=35
    ibn=bpy(:,2); ref=zeros(size(p,1),1); ref(ibn)=1;
    ref1=ref(t(:,1)); ref2=ref(t(:,2)); ref3=ref(t(:,3)); ref4=ref(t(:,4));
    rr=ref1+ref2+ref3+ref4;
    irr=find(rr==3); tn=t(irr,:); 
    tref=[ref1(irr) ref2(irr) ref3(irr) ref4(irr)]; [~,iref]=sort(tref,2);
    tb=zeros(size(tn));  
    for i=1:size(tn,1), tb(i,:)=tn(i,iref(i,:)); end
    tb(:,1)=[];  
    ibcn=union(tb(:,1),union(tb(:,2),tb(:,3)));
    clear bpx pby bpz ref1 ref2 ref3 ref4 rr irr tn tref iref tn
else
    error('The value of d is not valid')
end

np=size(p,1); nt=size(t,1);
fprintf('Dimension             :%3d\n',d)
fprintf('Number of vertices    :%5d\n',np)
fprintf('Number of elements    :%5d\n',nt)
fprintf('Number of d.o.f.      :%5d\n',d*np)

 
% Right-hand sides
f1=zeros(np,1);  
g1=zeros(np,1); g2=g1;
if (d==2)
   f2=fz*ones(np,1);
   g2(ibcn)=gz;
else
    f2=zeros(np,1); f3=fz*ones(np,1);
    g2=g1; g3=g1; g3(ibcn)=gz; 
end


% matrices & vecteurs
if (d==2)
    R=kelas2dstf(p,t,E,nu);
    f=kelas2drhs(p,t,f1,f2);
    g=kelas2drhsn(p,eneum,g1,g2);
else
    fprintf('Technique of Koko: ')
    tic; R=kelas3dstf(p,t,E,nu); toc;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% added by Jan Valdman 
%     add_paths; 
    
    fprintf('Technique of Rahman and Valdman: ')
    mu=E/(2*(1+nu));lambda=E*nu/((1+nu)*(1-2*nu)); 
    tic; R2=stiffness_matrixP1_3D_elasticity(t,p,lambda,mu); toc;
    %normest(R-R2)
    
    fprintf('Technique of Cermak, Sysala and Valdman: ')
    shear = E/(2*(1+nu)); % shear modulus
    bulk = E/(3*(1-2*nu)); % bulk modulus
    tic;  
    [Xi, WF] = quadrature_volume('P1');
    [HatP,DHatP1,DHatP2,DHatP3] = local_basis_volume('P1', Xi);
    n_e=size(t',2);     % number of elements
    n_q=length(WF);           % number of quadratic points
    n_int = n_e*n_q ;         % total number of integrations points
    shear =shear*ones(1,n_int); bulk=bulk*ones(1,n_int);
    [R3,WEIGHT]=elastic_stiffness_matrix(t',p',shear,bulk,DHatP1,DHatP2,DHatP3,WF); 
    %normest(R-R3)
    toc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    f=kelas3drhs(p,t,f1,f2,f3);
    g=kelas3drhsn(p,tb,g1,g2,g3);
end

% Dirichlet conditions (penalization)
R(ibcd,ibcd)=R(ibcd,ibcd)+10^15*speye(length(ibcd));
b=f+g; b(ibcd)=0;
 
% Solution (Gaussian elimination)
u=R\b;  

% vizualisation
if (d==2)
    s=kelas2dvmes(p,t,u,E,nu);
    kelas2dshow(p,t,u,100,s), axis off
else
    s=kelas3dvmes(p,t,u,E,nu);
    kelas3dshow(p,t,u,100,s), axis off
end
