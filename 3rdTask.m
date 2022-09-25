clear
clc
dt = 2.0;
k1 = 0.4;
k2 = 0.5;
k3 = 0.25;
k4 = 0.3;
for i=1:100
   C_A(i)=0.0;
   C_B(i)=0.0;
   C_C(i)=0.0;
   C_D(i)=0.0;
   C_E(i)=0.0;
   C_F(i)=0.0;
   C_G(i)=0.0;
   C_A_diff(i) = 0.0;
   C_B_diff(i) = 0.0;
   C_C_diff(i) = 0.0;
   C_D_diff(i) = 0.0;
   C_E_diff(i)=0.0;
   C_F_diff(i)=0.0;
   C_G_diff(i)=0.0;
end
C_A(1)=0.0625;
C_B(1)=0.138;
C_C(1)=0.15;
C_D(1)=0.0;
C_E(1)=0.0;
C_F(1)=0.0;
C_G(1)=0.0;
C_A_diff(1) = 0.0;
C_B_diff(1) = 0.0;
C_C_diff(1) = 0.0;
C_D_diff(1) = 0.0;
C_E_diff(1)=0.0;
C_F_diff(1)=0.0;
C_G_diff(1)=0.0;
for i = 1 : 99
   
   C_A(i+1) = C_A(i) + dt * (-k1 * C_A(i)*C_B(i)*(C_C(i)*C_C(i))+k3*C_E(i));
   C_B(i+1) = C_B(i) + dt * (-k1 * C_A(i)* C_B(i) *(C_C(i))^2 + k2* C_D(i)-2*k4*C_B(i)*C_B(i));
   C_C(i+1) = C_C(i) + dt * (-2*k1 * C_A(i)*C_B(i) *(C_C(i))^2 + k2*C_D(i)+ k3*C_E(i));
   C_D(i+1) = C_D(i) + dt * (k1*C_A(i) * C_B(i) * (C_C(i))^2 - k2*C_D(i));
   
    C_E(i+1) = C_E(i) + dt * (k2*C_D(i) - k3*C_E(i));
    C_F(i+1)=C_F(i)+dt*(2*k4*(C_B(i))^2);
    C_G(i+1)=C_G(i)+dt*(k4*(C_B(i))^2);
end
time = dt * (1:100);
fig = figure();
plot(time, C_A, 'b-.', 'LineWidth', 2)
hold on
plot(time, C_B, 'r', 'LineWidth', 2)
hold on
plot(time, C_C, 'y--', 'LineWidth', 2)
hold on
plot(time, C_D, 'g--', 'LineWidth', 2)
hold on
plot(time, C_E, 'b--.', 'LineWidth', 2)
hold on
plot(time, C_F, 'r--', 'LineWidth', 2)
hold on
plot(time, C_G, 'y-.', 'LineWidth', 2)
hold on
 
set(gca, 'FontSize', 18)
set(fig, 'color', 'white')
grid on
xlabel('time [sec]')
ylabel('Concentration [mole/cm^3]')
ylim([min(0.0) max(0.8)])
legend('C(A)','C(B)','C(C)','C(D)', 'C(E)', 'C(F)', 'C(G)');
 
%By Runge Kutta method:
 
dt = 2.0;
k1 = 0.4;
k2 = 0.5;
k3 = 0.25;
k4 = 0.3;
for i=1:100
   C_A2(i)=0.0;
   C_B2(i)=0.0;
   C_C2(i)=0.0;
   C_D2(i)=0.0;
   C_E2(i)=0.0;
   C_F2(i)=0.0;
   C_G2(i)=0.0;
end
C_A2(1)=0.0625;
C_B2(1)=0.138;
C_C2(1)=0.15;
C_D2(1)=0.0;
C_E2(1)=0.0;
C_F2(1)=0.0;
C_G2(1)=0.0;
for i = 1 : 99
% C_A
   k1_1 = -k1 * C_A2(i)*C_B2(i)*(C_C2(i))^2+k3*C_E2(i);
   k2_2 = -k1 * (C_A2(i) + dt * k1_1 / 2)*C_B2(i)*(C_C2(i))^2+k3*C_E2(i);
   k3_3 = -k1 * (C_A2(i) + dt * k2_2 / 2)*C_B2(i)*(C_C2(i))^2+k3*C_E2(i);
   k4_4 = -k1 * (C_A2(i) + dt * k3_3 / 2)*C_B2(i)*(C_C2(i))^2+k3*C_E2(i);
   deltaY1 = dt / 6 * (k1_1 + 2 * k2_2 + 2 * k3_3 + k4_4);
   C_A2(i+1) = C_A2(i) + deltaY1;
% C_B
   k1_1 = (-k1 * C_A2(i)* C_B2(i) *(C_C2(i))^2 + k2* C_D2(i)-2*k4*C_B2(i)*C_B2(i));
   k2_2 = (-k1 * C_A2(i)* (C_B2(i) + dt * k1_1 / 2) *(C_C2(i))^2 + k2* C_D2(i)-2*k4*(C_B2(i) + dt * k1_1 / 2)*(C_B2(i) + dt * k1_1 / 2));
   k3_3 = (-k1 * C_A2(i)* (C_B2(i) + dt * k2_2 / 2) *(C_C2(i))^2 + k2* C_D2(i)-2*k4*(C_B2(i) + dt * k2_2 / 2)*(C_B2(i) + dt * k2_2 / 2));
   k4_4 = (-k1 * C_A2(i)* (C_B2(i) + dt * k3_3 / 2) *(C_C2(i))^2 + k2* C_D2(i)-2*k4*(C_B2(i) + dt * k3_3 / 2)*(C_B2(i) + dt * k3_3 / 2));
   deltaY2 = dt / 6 * (k1_1 + 2 * k2_2 + 2 * k3_3  + k4_4);
   C_B2(i+1) = C_B2(i) + deltaY2;
% C_C
   k1_1 =  -2*k1 * C_A2(i)*C_B2(i) *(C_C2(i))^2 + k2*C_D2(i)+ k3*C_E2(i);
   k2_2 =  -2*k1 * C_A2(i)*C_B2(i) *((C_C2(i) + dt * k1_1 / 2))^2 + k2*C_D2(i)+ k3*C_E2(i);
   k3_3 =  -2*k1 * C_A2(i)*C_B2(i) *((C_C2(i) + dt * k2_2 / 2))^2 + k2*C_D2(i)+ k3*C_E2(i);
   k4_4 =  -2*k1 * C_A2(i)*C_B2(i) *((C_C2(i) + dt * k3_3 / 2))^2 + k2*C_D2(i)+ k3*C_E2(i);
   deltaY3 = dt / 6 * (k1_1 + 2 * k2_2 + 2 * k3_3  + k4_4);
   C_C2(i+1) = C_C2(i) + deltaY3;
% C_D
   k1_1 = k1*C_A2(i) * C_B2(i) * (C_C2(i))^2 - k2*C_D2(i);
   k2_2 = k1*C_A2(i) * C_B2(i) * (C_C2(i))^2 - k2*(C_D2(i) + dt * k1_1 / 2);
   k3_3 = k1*C_A2(i) * C_B2(i) * (C_C2(i))^2 - k2* (C_D2(i) + dt * k2_2 / 2);
   k4_4 = k1*C_A2(i) * C_B2(i) * (C_C2(i))^2 - k2* (C_D2(i) + dt * k3_3 / 2);
   deltaY4 = dt / 6 * (k1_1 + 2 * k2_2 + 2 * k3_3  + k4_4);
   C_D2(i+1) = C_D2(i) + deltaY4;
% C_E
   k1_1 = k2*C_D2(i) - k3*C_E2(i);
   k2_2 = k2*C_D2(i) - k3*(C_E2(i) + dt * k1_1 / 2);
   k3_3 = k2*C_D2(i) - k3*(C_E2(i) + dt * k2_2 / 2);
   k4_4 = k2*C_D2(i) - k3*(C_E2(i) + dt * k3_3 / 2);
   deltaY5 = dt / 6 * (k1_1 + 2 * k2_2 + 2 * k3_3  + k4_4);
   C_E2(i+1) = C_E2(i) + deltaY5;
   
 % C_F
   k1_1 =(2*k4*(C_B2(i))^2);
   k2_2 = (2*k4*(C_B2(i))^2);
   k3_3 = (2*k4*(C_B2(i))^2);
   k4_4 = (2*k4*(C_B2(i))^2);
   deltaY6 = dt / 6 * (k1_1 + 2 * k2_2 + 2 * k3_3  + k4_4);
   C_F2(i+1) = C_F2(i) + deltaY6;
  %C_G
   k1_1 =(k4*(C_B2(i))^2);
   k2_2 = (k4*(C_B2(i))^2);
   k3_3 = (k4*(C_B2(i))^2);
   k4_4 = (k4*(C_B2(i))^2);
   deltaY7 = dt / 6 * (k1_1 + 2 * k2_2 + 2 * k3_3  + k4_4);
   C_G2(i+1) = C_G2(i) + deltaY7;
     
   C_A_diff(i) = abs(C_A(i) - C_A2(i));
   C_B_diff(i) = abs(C_B(i) - C_B2(i));
   C_C_diff(i) = abs(C_C(i) - C_C2(i));
   
   C_D_diff(i) = abs(C_D(i) - C_D2(i));
   C_E_diff(i) = abs(C_E(i) - C_E2(i));
   C_F_diff(i) = abs(C_F(i) - C_F2(i));
   C_G_diff(i) = abs(C_G(i) - C_G2(i));
   fprintf('%.0f: %.7f - %.7f = %.7f\n', i, C_A(i), C_A2(i), abs(C_A(i) - C_A2(i)));
end
time = dt * (1:100);
fig = figure();
set(fig, 'color', 'white')
plot(time, C_A2, 'b-.', 'LineWidth', 2)
hold on
plot(time, C_B2, 'r', 'LineWidth', 2)
hold on
plot(time, C_C2, 'y--', 'LineWidth', 2)
hold on
plot(time, C_D2, 'g--', 'LineWidth', 2)
hold on
plot(time, C_E2, 'b-', 'LineWidth', 2)
hold on
plot(time, C_F2, 'r--', 'LineWidth', 2)
hold on
plot(time, C_G2, 'y-.', 'LineWidth', 2)
hold on
set(gca, 'FontSize', 18)
set(fig, 'color', 'white')
grid on
xlabel('time [sec]')
ylabel('Concentration [mole/cm^3]')
ylim([min(0.0) max(0.8)])
legend('C(A)','C(B)','C(C)', 'C(D)', 'C(E)', 'C(F)', 'C(G)');
 
fig = figure();
set(fig, 'color', 'white')
plot(time, C_A_diff, 'b-.', 'LineWidth', 2)
hold on
plot(time, C_B_diff, 'k--', 'LineWidth', 2)
hold on
plot(time, C_C_diff, 'y', 'LineWidth', 2)
hold on
plot(time, C_D_diff, 'g', 'LineWidth', 2)
hold on
plot(time, C_E_diff, 'b--', 'LineWidth', 2)
hold on
plot(time, C_F_diff, 'c-.', 'LineWidth', 2)
hold on
plot(time, C_G_diff, 'y--', 'LineWidth', 2)
hold on
set(gca, 'FontSize', 18)
set(fig, 'color', 'white')
ax = gca;
ax.YAxis.Exponent = 0;
%ylim([0 0.25])
grid on
xlabel('time [sec]')
ylabel('Concentration [mole/cm^3]')
ylim([min(0.0) max(0.8)])
legend('C(A)','C(B)','C(C)', 'C(D)','C(E)', 'C(F)', 'C(G)');
