% Simulating voltages on surface of material. left and right at 1V and top
% and bottom at 0V.

clear all
clearvars
clearvars -GLOBAL
close all
format shorte

numIt = 10; % number of iterations to take
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

            if (i == 1) || (i == nx) % left and right sides = 1V

                G(n,n) = 1;
                B(n) = 1;

            elseif (j == ny) || (j == 1) % right side = 0V

                G(n,n) = 1;

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
    
    surf(VmapInv) % plot surface
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

subplot(2,1,2)
quiver(Ex,Ey)

Jx = cond*Ex;
Jy = cond*Ey;

% calculate total current
Cinit = sum(Jx(1,:));
Cfin = sum(Jx(nx,:));
totCurr = (Cinit + Cfin) / 2.0;