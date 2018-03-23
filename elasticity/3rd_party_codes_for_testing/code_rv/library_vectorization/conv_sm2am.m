function ama = conv_sm2am (smx,nz)
% ama: ama(1:ny,1:nx,1:nz)
% smx: smx(1:nk,1:nx)
% amb: amb(1:nk,1:ny,1:nz)

[nk,nx]    = size(smx);

ama=zeros(nk,nx,nz);

for row = 1:nk
    for col = 1:nx
        ama(row,col,:)=smx(row,col); 
    end
end


