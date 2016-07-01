% 3Dtrial.m

clear
clc
close('all');

Num_ts = 1000;
ts_freq = 10;
plot_freq = 50;

dynamics = 3;
% 1 = LBGK
% 2 = TRT 
% 3 = MRT <--- D3Q15 and D3Q19 only

Lx = 4;  %Length along channel
Ly = 1;  %Width of channel
Lz = 1; %Height of channel

x_c = 2.0; %Location of pile centroid
y_c = 0.5; %Center of channel
r_c = 0.1; %radius of pile
z_b = 0; %Bottom of pile
h_c = Lz; %Pile height equals channel height
z_t = z_b + h_c;

%Scour Pit Dimensions
ellip_a = 4*2*r_c; %4x diameter of pile
ellip_b = 2*2*r_c; %2x diam of pile
ellip_c = 1*2*r_c; %1x diam of pile
ellip_x = x_c + r_c; %center located at back of pile in x-dir
ellip_y = y_c; %center located in center of pile
ellip_z = ellip_c; %bottom of ellipse at zm

rho_p = 1260; %kg/m^3 (glycerin)
nu_p = 1.49/rho_p; %m^2/s

u_in = 0.44; %m/s velocity into the channel

Re = u_in*2*r_c/nu_p; %Using existing 

%Reynolds number (open channel) = 4*hyd_rad*u_in/nu
%hyd_rad = Lz; %Large open channels - r_h simplifies to A/width = depth
%Re = 4*hyd_rad*u_in/nu_p; 
fprintf('Reynolds Number = %g \n',Re);

%Dimensionless units
Lo = 2*r_c; %Characteristic length
Uo = u_in;
To = Lo/Uo;

Ld =1; 
Td = 1;
Ud = 1;
nu_d = 1/Re;

%convert to LBM units
dt = 0.0025;
Nx_divs = 11; dx = 1/(Nx_divs - 1);
Ny_divs = 11; dy = 1/(Ny_divs - 1);
Nz_divs = 11; dz = 1/(Nz_divs - 1);
u_lbm = (dt/dx)*Ud;
nu_lbm = (dt/(dx^2))*nu_d;
omega = 1/(3*nu_lbm + 0.5);

u_conv_fact = (dt/dx)*(To/Lo);
t_conv_fact = (dt*To);
l_conv_fact = (dt*Lo);

rho_lbm = rho_p;
rho_out = rho_lbm;

%make lattice
xm = 0; xp = Lx;
ym = 0; yp = Ly;
zm = 0; zp = Lz;

Nx = ceil((Nx_divs-1)*(Lx/Lo))+1;
Ny = ceil((Ny_divs-1)*(Ly/Lo))+1;
Nz = ceil((Nz_divs-1)*(Lz/Lo))+1;

x_space = linspace(xm,xp,Nx);
y_space = linspace(ym,yp,Ny);
z_space = linspace(zm,zp,Nz);

[X,Y,Z] = meshgrid(x_space,y_space,z_space);

ind_init = 1:(Nx*Ny*Nz); ind = reshape(ind_init,[Ny Nx Nz]);

gridcoord = [X(:) Y(:) Z(:)];

% elements will be constructed using the same node-ordering in x-y plane as
% with 2D elements, with the 5-6-7-8 nodes occuring on the next higher z
% level

ind1 = ind(1:(end-1),1:(end-1),1:(end-1));
ind1 = ind1(:);
ind2 = ind(1:(end-1),2:end,1:(end-1)); ind2=ind2(:);
ind3 = ind(2:end,2:end,1:(end-1)); ind3=ind3(:);
ind4 = ind(2:end,1:(end-1),1:(end-1)); ind4=ind4(:);
ind5 = ind(1:(end-1),1:(end-1),2:end); ind5=ind5(:);
ind6 = ind(1:(end-1),2:end,2:end); ind6=ind6(:);
ind7 = ind(2:end,2:end,2:end); ind7=ind7(:);
ind8 = ind(2:end,1:(end-1),2:end); ind8=ind8(:);

nodes=[ind1,ind2,ind3,ind4,ind5,ind6,ind7,ind8];

faces.zy_m = find(gridcoord(:,1)==xm);
faces.zy_p = find(gridcoord(:,1)==xp);
faces.zx_m = find(gridcoord(:,2)==ym);
faces.zx_p = find(gridcoord(:,2)==yp);
faces.xy_m = find(gridcoord(:,3)==(zm:0.5*ellip_c));
faces.xy_p = find(gridcoord(:,3)==zp);

fprintf('dx = %g, dy = %g, dz = %g \n',x_space(2) - x_space(1),...
    y_space(2)-y_space(1), z_space(2)-z_space(1));

[w,ex,ey,ez,bb_spd]=D3Q15_lattice_parameters();

streamTgtMat = genStreamTgtVec3D(Nx,Ny,Nz,ex,ey,ez);

nnodes = Nx*Ny*Nz;

fprintf('Number of Nodes = %g \n',nnodes);

%D3Q15
numSpd = 15;

M = getMomentMatrix('D3Q15');

S = omega.*eye(numSpd);

omega_op = M\(S*M);

%Cylinder Dimensions - Circular cylinder with variable height

cyl_list = find(((gridcoord(:,1)-x_c).^2 + (gridcoord(:,2)-y_c).^2 ...
    < (r_c)^2) & z_b<=gridcoord(:,3)<=z_t);

%Solid Nodes - open channel, thick bottom allows for scour pit development.
solidnodes = find(gridcoord(:,3)<=ellip_c);

%This is the ellipsoid that makes up the scour pit.
scourpit = find(((gridcoord(:,1) - ellip_x).^2/ellip_a.^2) + ((gridcoord(:,2)-y_c).^2/ellip_b.^2)...
    +((gridcoord(:,3) - ellip_z).^2/ellip_c.^2) <= 1);

%Remove ellipsoid nodes from bottom of channel.
solidnodes = setxor(solidnodes,intersect(solidnodes,scourpit));

%solidnodes = unique(solidnodes); %Needed if solid boundaries
%located on sides or top.  Save this.

% add the cylindrical obstacle to the solid node list
solidnodes = [solidnodes; cyl_list];

inletnodes = faces.zy_m; 
inletnodes = setxor(inletnodes,intersect(inletnodes,solidnodes)); % eliminate solid nodes from inlet nodes
outletnodes = faces.zy_p;
outletnodes = setxor(outletnodes,intersect(outletnodes,solidnodes)); %eliminate solid nodes from outletnodes

% Potentially use later to identify new boundary for full-slip.
% topnodes = faces.xy_p;
% topnodes = setxor(topnodes,intersect(topnodes,solidnodes));

%Set Inlet Velocities.  Constant velocity along inlet/outlet.
ux_in = u_lbm;
uy_in = zeros(length(ux_in),1);
uz_in = zeros(length(ux_in),1);

ux_out = ux_in;
uy_out = uy_in;
uz_out = uz_in;

%Visualization Settings
% tag some nodes for visualization
% xy-plane along the centerline to get an idea of the overall flow pattern
z_pln = z_space(ceil(Nz/2));
vis_nodes = find(gridcoord(:,3)==z_pln);

% zy-plane at mid-channel to get a view of the flow profile
zy_pln = ind(:,ceil(Nx/2),:);
zy_pln = reshape(zy_pln,[Ny Nz]);

y_pln = y_space(ceil(Ny/2));
vis_nodes_y = find(gridcoord(:,2)==y_pln);

%Initialization
[fIn,fOut,rho,ux,uy,uz]=Initialize_F_zero3D(gridcoord,ex,ey,ez,w,rho_lbm);

switch dynamics
    
    case 1% BGK
        fEq = zeros(nnodes,numSpd);
    case 2 % TRT
        fEq = zeros(nnodes,numSpd);
        fNEq = zeros(nnodes,numSpd);
        fEven = zeros(nnodes,numSpd);
        fOdd = zeros(nnodes,numSpd);
        
    case 3 % MRT
        fEq = zeros(nnodes,numSpd);
        M = getMomentMatrix('D3Q15');
        S = getEwMatrixMRT('D3Q15',omega);
        omega_op = M\(S*M);
end

fprintf('Number of Lattice-points = %d.\n',nnodes);
fprintf('Number of time-steps = %d. \n',Num_ts);
%fprintf('Predicted execution time = %g.\n', predicted_ex_time);

fprintf('LBM viscosity = %g. \n',nu_lbm);
fprintf('LBM relaxation parameter (omega) = %g. \n',omega);
fprintf('LBM flow Mach number = %g. \n',u_lbm);

input_string = sprintf('Do you wish to continue? [Y/n] \n');

run_dec = input(input_string,'s');

if ((run_dec ~= 'n') && (run_dec ~= 'N'))
    
    fprintf('Here we go... \n');
    
  % commence time stepping
    
    tic;
    for ts = 1:Num_ts
        
        % say something comforting about my progress...
        if(mod(ts,ts_freq)==0)
            fprintf('Executing time step number %d.\n',ts);
        end
        
        % compute density
        rho = sum(fIn,2);
        
        % compute velocities
        ux = (fIn*ex')./rho;
        uy = (fIn*ey')./rho;
        uz = (fIn*ez')./rho;
        
        % set macroscopic and Microscopic Dirichlet-type boundary
        % conditions
        
        % macroscopic BCs
        ux(inletnodes)=ux_in;
        uy(inletnodes)=uy_in;
        uz(inletnodes)=uz_in;
        
        ux(outletnodes)=ux_out;
        uy(outletnodes)=uy_out;
        uz(outletnodes)=uz_out;
        
               
        
        % microscopic BCs
        fIn(inletnodes,:)=velocityBC_3D(fIn(inletnodes,:),w,ex,ey,ez,...
            ux_in,uy_in,uz_in);
        fIn(outletnodes,:)=velocityBC_3D(fIn(outletnodes,:),w,ex,ey,ez,...
            ux_out,uy_out,uz_out);
        % Collide
        switch dynamics
            
            case 1 % LBGK
                for i = 1:numSpd
                    cu = 3*(ex(i)*ux+ey(i)*uy+ez(i)*uz);
                    fEq(:,i)=w(i)*rho.*(1+cu+(1/2)*(cu.*cu) - ...
                        (3/2)*(ux.^2 + uy.^2+uz.^2 ));
                    fOut(:,i)=fIn(:,i)-omega*(fIn(:,i)-fEq(:,i));
                end
                
            case 2 % TRT
                
                % find even part and odd part of off-equilibrium
                % for all speeds, then relax those parts...
                
                % compute fEq
                for i = 1:numSpd
                    cu = 3*(ex(i)*ux+ey(i)*uy+ez(i)*uz);
                    fEq(:,i)=w(i)*rho.*(1+cu+(1/2)*(cu.*cu) - ...
                        (3/2)*(ux.^2 + uy.^2+uz.^2 ));
                end
                
                % compute fNEq
                fNEq = fEq - fIn;
                
                % compute odd and even part of fNEq
                for i= 1:numSpd
                    fEven(:,i) = (1/2)*(fNEq(:,i)+fNEq(:,bb_spd(i)));
                    fOdd(:,i) = (1/2)*(fNEq(:,i)-fNEq(:,bb_spd(i)));
                    
                end
                
                % relax on these even and odd parts
                fOut = fIn + omega*fEven + fOdd;
                % omega for odd part is equal to 1.
                
            case 3 % MRT
                % compute fEq
                for i = 1:numSpd
                    cu = 3*(ex(i)*ux+ey(i)*uy+ez(i)*uz);
                    fEq(:,i)=w(i)*rho.*(1+cu+(1/2)*(cu.*cu) - ...
                        (3/2)*(ux.^2 + uy.^2+uz.^2 ));
                end
                % collide
                fOut = fIn - (fIn - fEq)*omega_op;
                
                
        end
        
         % bounce-back
        for i = 1:numSpd
            fOut(solidnodes,i)=fIn(solidnodes,bb_spd(i));
        end
        
        
        % stream
        for i = 1:numSpd
            fIn(streamTgtMat(:,i),i)=fOut(:,i);
        
       
        end
        
        % visualization        
        if(mod(ts,plot_freq)==0)
           
             figure(1)
            ux_vp = ux(vis_nodes);
            uy_vp = uy(vis_nodes);
            uz_vp = uz(vis_nodes);
            u_vp = sqrt(ux_vp.^2+uy_vp.^2+uz_vp.^2)./u_conv_fact;
            u_vp =reshape(u_vp,[Ny Nx]);
            imagesc(u_vp);
            colorbar
            title('Midplane velocity profile')
            axis equal off
            drawnow
            
            figure(2)
            ux_v_2 = ux(vis_nodes_y);
            ux_v_2 =reshape(ux_v_2,[Nx Nz]);
            imagesc(ux_v_2);
            colorbar
            title('Longitudinal View - u_x only')
            axis equal 
            xlabel('z')
            ylabel('x')
            camroll(90)
            drawnow
            
         figure(3)
            uz_v_2 = uz(vis_nodes_y);
            uz_v_2 =reshape(uz_v_2,[Nx Nz]);
            imagesc(uz_v_2);
            colorbar
            title('Longitudinal View - u_z only')
            axis equal 
            xlabel('z')
            ylabel('x')
            camroll(90)
            drawnow
            
            
            figure(4)
            ux_cp = ux(zy_pln(:));
            uy_cp = uy(zy_pln(:));
            uz_cp = uz(zy_pln(:));
            u_cp = sqrt(ux_cp.^2+uy_cp.^2+uz_cp.^2)./u_conv_fact;
            u_cp = reshape(u_cp,[Ny Nz]);
            mesh(u_cp)
            view(3)
            title('Mid-channel velocity profile')
            drawnow
            
          
            
            
        end
        
    end
    ex_time = toc;
    
    fprintf('Lattice-point updates per second = %g.\n',...
        Num_ts*nnodes/ex_time);
    
    
    
    
else
    fprintf('Run aborted.  Better luck next time!\n');
end