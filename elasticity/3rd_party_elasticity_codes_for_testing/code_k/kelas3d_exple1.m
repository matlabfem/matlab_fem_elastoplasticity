%
% Linear elasticity  2D/3D
% 
% 2D: \Omega=(-5 5)*(0 35)
%     f=(0,-.75)
%     g=(0,10)
% 3D: \Omega=(-5 5)*(0 35)*(-5 5)
%     f=(0,0,-.75)
%     g=(0,0,10)
%----------------------------------------------------------
E=30000; nu=0.4;
d=2;
nx=11;
ny=36;
fx=0; fy=0; fz=-.75;
gx=0; gy=0; gz=10;

% mesh generation
if (d==2)
    ax=-5; bx=5; ay=0; by=35;  
    [p,t,bpx,bpy]=kpde2dumsh(ax,bx,ay,by,nx,ny);
    ibcd=[2*bpy(:,1)-1; 2*bpy(:,1)-1];
    ibn=bpy(:,2); x=p(ibn,1); [~,ix]=sort(x); ibn=ibn(ix); ne=length(ibn);
    eneum=[ibn(1:ne-1) ibn(2:ne)];
    ibcn=ibn;
    clear bpx bpy ibn x ix
elseif (d==3)
    ax=-5; bx=5; ay=0; by=35; az=-5; bz=5;
    nz=nx;
    [p,t,bpx,bpy,bpz]=kpde3dumsh(ax,bx,ay,by,az,bz,nx,ny,nz);
    ibcd=[3*bpy(:,1)-2; 3*bpy(:,1)-1; 3*bpy(:,1)];
    % Neumann condition at y=35
    ibn=bpy(:,2); ref=zeros(np,1); ref(ibn)=1;
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
   f2=-fz*ones(np,1);
   g2(ibcn)=gz;
else
    f2=zeros(np,1); f3=-fz*ones(np,1);
    g2=g1; g3=g1; g3(ibcn)=gz; 
end


% matrices & vecteurs
if (d==2)
    R=kelas2d(p,t,E,nu);
    f=kelas2drhs(p,t,f1,f2,f3);
    g=kelas2drhsn(p,eneum,g1,g2);
else
    R=kelas3d(p,t,E,nu);
    f=kelas3drhs(p,t,f1,f2,f3,);
    g=kelas3drhsn(p,tb,g1,g2);
end

% Dirichlet conditions (penalization)
R(ibcd,ibcd)=A(ibcd,ibcd)+10^20*speye(length(ibcd));
b=f+g; b(ibcd)=0;
 
% Solution (Gaussian elimination)
ui=A\b; u=zeros(nn,1); u(inodes)=ui;

% vizualisation
if (d==2)
    s=kelas2dvmes(p,t,u,E,nu);
    kelas2dshow(p,t,100,s)
else
    s=kelas3dvmes(p,t,u,E,nu);
    kelas3dshow(p,t,u,100,s)
end
