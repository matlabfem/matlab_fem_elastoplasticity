function av=conv_ma2av(ma)
%converts matrix to array of vectors
% ma: ma(1:nx, 1:ny)
% av: av(1:y,1,1:nx)

 [nx,ny] = size(ma);
  av=zeros(ny,1,nx);
  
  for i=1:ny
      av(i,1,:)=ma(:,i);
  end

end
