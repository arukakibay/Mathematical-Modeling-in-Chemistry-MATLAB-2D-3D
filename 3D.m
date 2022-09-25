     clc;
clear;
x_min = 0;
x_max = 15;
 
N = x_max * 10;
M = x_max * 10;
 
dx = (x_max - x_min)/N; 
dy = (x_max - x_min)/N;
 
figure('units','normalized','outerposition',[0 0 1 1]);
 
x = linspace(0,x_max,N);
y = linspace(0,x_max,M);
[X,Y] = meshgrid(x,y);
T_0 = 0;
 
A_old = ones(N,M)*T_0;
A_new = ones(N,M)*T_0;
B_old = ones(N,M)*T_0;
B_new = ones(N,M)*T_0;
C_old = ones(N,M)*T_0;
C_new = ones(N,M)*T_0;
D_old = ones(N,M)*T_0;
D_new = ones(N,M)*T_0;
 
for i = 1:N
    for j = 1:M
        if sqrt((7.5 - i * dx)^2 + (7.5 - j * dy)^2) <= 1 
            A_old(i,j) = 7;
        end
    end
end
 
 
 
t = 0; t_end = 24;
dt = 2 * dx * dy;
 
a_1 = 5 * dx * dx; 
a_2 = 5 * dx * dx; 
b_1 = 0.5; 
b_2 = 0.5; 
b_3 = 4 * dx * dx; 
b_4 = 4 * dx * dx; 
c_1 = 3 * dx * dx; 
c_2 = 3 * dx * dx; 
d_1 = 1; 
d_2 = 1; 
 
k = 0.1;
 
while t < t_end   
    t = t + dt;
    for i = 2:N-1
        for j = 2:M-1
            A_new(i,j) = A_old(i,j) + dt * (-2 * k * A_old(i,j)^2 * B_old(i,j)                                                                                 + a_1 / (dx * dx) * (A_old(i+1,j) - 2 * A_old(i,j) + A_old(i-1,j)) + a_2 / (dy * dy) * (A_old(i,j+1) - 2 * A_old(i,j) + A_old(i,j-1)) );
            B_new(i,j) = B_old(i,j) + dt * (    -k * A_old(i,j)^2 * B_old(i,j) - b_1 / dx * (B_old(i,j) - B_old(i-1,j)) - b_2 / dx * (B_old(i,j) - B_old(i,j-1)) + b_3 / (dx * dx) * (B_old(i+1,j) - 2 * B_old(i,j) + B_old(i-1,j)) + b_4 / (dy * dy) * (B_old(i,j+1) - 2 * B_old(i,j) + B_old(i,j-1)) );
            C_new(i,j) = C_old(i,j) + dt * (     k * A_old(i,j)^2 * B_old(i,j)                                                                                 + c_1 / (dx * dx) * (C_old(i+1,j) - 2 * C_old(i,j) + C_old(i-1,j)) + c_2 / (dy * dy) * (C_old(i,j+1) - 2 * C_old(i,j) + C_old(i,j-1)) );
            D_new(i,j) = D_old(i,j) + dt * ( 3 * k * A_old(i,j)^2 * B_old(i,j) - d_1 / dx * (D_old(i,j) - D_old(i-1,j)) - d_2 / dx * (D_old(i,j) - D_old(i,j-1))                                                                                                                                     );
        end
    end
    if t >= 0.5 && t <= 1.0
        for i = 1:N
            for j = 1:M
                if sqrt((7.5/2 - i * dx/2)^2 + (7.5/2 - j * dy/2)^2) <= 1 
                    B_new(i,j) = (t-0.5) * 3 * 2;
                end
            end
        end
    end
    A_old = A_new;
    B_old = B_new; 
    C_old = C_new; 
    D_old = D_new; 
    
    surf(X,Y,A_old);
    shading interp;
    hold on;    
    surf(X,Y,B_old);
    shading interp;
    hold on;    
    surf(X,Y,C_old);
    shading interp;
    hold on;    
    surf(X,Y,D_old);
    shading interp;
    hold off;
    %view(2);
    colormap(jet(90)), colorbar
    caxis([0 7])
    title({sprintf('\n time = %1.3f [sec]',t)},'FontSize',14)
    zlim([0 7])
    xlabel('X (cm)','FontSize',12);
    ylabel('Y (cm)','FontSize',12);
    zlabel('Concentration (mole/cm^{2})','FontSize',12);    
    axis([0 x_max 0 x_max 0 7]);
    drawnow;
end

