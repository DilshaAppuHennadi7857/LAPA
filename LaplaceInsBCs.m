% Simulating voltages on surface of material. Left side at 1V, right side
% at 0V, top and bottom both insultated.

clear all
clearvars
clearvars -GLOBAL
close all
format shorte

numIt = 8; % number of iterations to take
nx = 60;
ny = 60;
cond = 10; % approx. conductivity of silicon, assume uniform conductivity

% Create G and B vectors
G = zeros(nx*ny);
B = zeros(1,nx*ny);

% Run iterations and calculate voltages at each node
% Currently have a set number of iterative loops, but this can be adjusted
% to check for convergence
for iteration = 1:numIt
    for i = 1:nx
        for j = 1:ny

            n = j + (i-1)*ny; % node mapping

            if i == 1 % left side = 1V

                G(n,n) = 1;
                B(n) = 1;

            elseif i == nx % right side = 0V

                G(n,n) = 1;

            elseif j == ny % top side = insulated

                % only three resistors:
                % n(x-1,y), n(x+1,y), n(x,y-1)

                nxm = j + (i-2)*ny;
                nxp = j + i*ny;
                nym = (j-1) + (i-1)*ny;

                % assume uniform conductivity throughout material
                rxm = cond/2.0;
                rxp = cond/2.0;
                rym = cond/2.0;

                G(n,n) = -(rxm + rxp + rym);
                G(n,nxm) = rxm;
                G(n,nxp) = rxp;
                G(n,nym) = rym;

            elseif j == 1 % bottom side = insulated

                % only three resistors:
                % n(x-1,y), n(x+1,y), n(x,y+1)

                nxm = j + (i-2)*ny;
                nxp = j + i*ny;
                nyp = (j+1) + (i-1)*ny;

                % assume uniform conductivity throughout material
                rxm = cond/2.0;
                rxp = cond/2.0;
                ryp = cond/2.0;

                G(n,n) = -(rxm + rxp + ryp);
                G(n,nxm) = rxm;
                G(n,nxp) = rxp;
                G(n,nyp) = ryp;

            else % middle node

                nxm = j + (i-2)*ny;
                nxp = j + i*ny;
                nym = (j-1) + (i-1)*ny;
                nyp = (j+1) + (i-1)*ny;

                rxm = cond/2.0;
                rxp = cond/2.0;
                rym = cond/2.0;
                ryp = cond/2.0;

                G(n,n) = -(rxm + rxp + rym + ryp);
                G(n,nxm) = rxm;
                G(n,nxp) = rxp;
                G(n,nym) = rym;
                G(n,nyp) = ryp;

            end

        end
    end
    V = G\B';

    % Map voltages back into a matrix

    Vmap = zeros(nx,ny); % initialize matrix
    n = 0; % clear/reset node index n

    for i = 1:nx
        for j = 1:ny

            n = j + (i-1)*ny;
            Vmap(i,j) = V(n);

        end
    end

    % Mapping needs to be inverted to set the sides of the matrix to match the
    % appropriate BCs
    VmapInv = Vmap';
    
    figure(1)
    surf(VmapInv) % plot surface
    view(45,45) % adjust camera angle for better view
    pause(0.001);
end

% E = grad(V)
% note: this version doesn't seem to recognize grad()
Ex = [];
Ey = [];

for i = 1:nx
    for j = 1:ny
        
        % Calculate Ex
        
        if i == 1
            
            Ex(i,j) = VmapInv(i+1,j) - VmapInv(i,j);
            
        elseif i == nx
            
            Ex(i,j) = VmapInv(i,j) - VmapInv(i-1,j);
            
        else
            
            Ex(i,j) = (VmapInv(i+1,j) - VmapInv(i-1,j)) / 2.0;
            
        end
        
        % Calculate Ey
        
        if j == 1
            
            Ey(i,j) = VmapInv(i,j+1) - VmapInv(i,j);
            
        elseif j == ny
            
            Ey(i,j) = VmapInv(i,j) - VmapInv(i,j-1);
            
        else
            
            Ey(i,j) = (VmapInv(i,j+1) - VmapInv(i,j-1)) / 2.0;
            
        end
        
    end
end

Ex = -Ex;
Ey = -Ey;

figure(2)
subplot(2,1,1)
surf(Ex,Ey)
view(45,45) % adjust camera angle to match figure 1's view

subplot(2,1,2)
quiver(Ex,Ey)

Jx = cond*Ex;
Jy = cond*Ey;

% calculate total current
Cinit = sum(Jx(1,:));
Cfin = sum(Jx(nx,:));
totCurr = (Cinit + Cfin) / 2.0;

% figure(3)
% imboxfilt(VmapInv,3);
%
% Cannot get imboxfilt() to work
