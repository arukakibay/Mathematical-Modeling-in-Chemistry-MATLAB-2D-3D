clear
figure('units','normalized','outerposition',[0 0 1 1])
 
x_min = 0;
x_max = 15;
 
N = x_max * 10;
 
dx = (x_max - x_min)/N; % 0.1
x = x_min : dx : x_max;
 
t = 0;
t_max = 4;
dt = dx * dx;
n_steps = t_max/dt;
 
 
a = 5 * dx * dx;
b1 = 0.5;
b2 = 4 * dx * dx;
c = 3 * dx * dx;
d = 1;
 
k = 0.5;
 
for i = 1 : N+1
    A0(i) = 0;
    B0(i) = 0;
    C0(i) = 0;
    D0(i) = 0;
end
 
for i = 1: N+1     
    if dx * i <= 6.086 || dx * i >= 8.914
        A0(i) = 0 * dx * i;
    else
        A0(i) = sqrt(20-(dx * i-7.5)^2/0.1);
    end
end
 
 
A_old = A0;
A_new = A0;
B_old = B0;
B_new = B0;
C_old = C0;
C_new = C0;
D_old = D0;
D_new = D0;
 
 
for i = 1 : n_steps
    A_old(1)=A_old(2);
    A_old(N+1)=A_old(N);
    B_old(1)=B_old(2);
    B_old(N+1)=B_old(N);
    C_old(1)=C_old(2);
    C_old(N+1)=C_old(N);
    D_old(1)=D_old(2);
    D_old(N+1)=D_old(N);
    
    for i=2: N
        A_new(i) = A_old(i) + dt * (-2 * k * A_old(i)^2 * B_old(i)                                      + a  / (dx * dx) * (A_old(i+1) - 2 * A_old(i) + A_old(i-1)) );

        B_new(i) = B_old(i) + dt * (    -k * A_old(i)^2 * B_old(i) - b1 /  dx * (B_old(i) - B_old(i-1)) + b2 / (dx * dx) * (B_old(i+1) - 2 * B_old(i) + B_old(i-1)) );

        C_new(i) = C_old(i) + dt * (     k * A_old(i)^2 * B_old(i)                                      + c  / (dx * dx) * (C_old(i+1) - 2 * C_old(i) + C_old(i-1)) );

        D_new(i) = D_old(i) + dt * ( 3 * k * A_old(i)^2 * B_old(i) - d  /  dx * (D_old(i) - D_old(i-1))                                                             );
    end
    
    if t >= 0.5 && t <= 1.0
        B_new = (t-0.5)*4*exp(-3*(x - 7.5).^2);
    end
    
    t=t+dt;
    A_old=A_new;
    B_old=B_new;
    C_old=C_new;
    D_old=D_new;
    
    plot(x,A_old,'b-','LineWidth',3);
    hold on
    plot(x,B_old,'r-','LineWidth',3);
    hold on
    plot(x,C_old,'g-','LineWidth',3);
    hold on
    plot(x,D_old,'y--','LineWidth',3);
    hold off
    axis([x_min x_max -0.5 10])
    xlabel('X','FontSize',16)
    ylabel('U(x,t)','FontSize',16)
    title(sprintf('time =%1.3f',t),'FontSize',16)
    legend('A','B','C', 'D');
     
    pause(dt);
end
