function stv = genStreamTgtVec3D(Nx,Ny,Nz,ex,ey,ez)

ind = 1:(Nx*Ny*Nz); ind = reshape(ind,[Ny Nx Nz]);

numSpd=length(ex);
nnodes = Nx*Ny*Nz;

stv = zeros(nnodes,numSpd);

for spd = 1:numSpd
   t = circshift(ind,[-ey(spd) -ex(spd) -ez(spd)]); t = t(:);
   stv(:,spd)=t;
     
end

    %ind_x = ind(2); ind_y = ind(1); ind_z = ind(3);
        %ind_z(ind_z==Nz)=0;
        
       % ind = [ind_x ind_y ind_z];