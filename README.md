% Preparación
clear all
clc
close all

% Información Conocida
global Lab
global Lbc
global Lcd
global Lda
global cranck_angle

Lab = 4183.417; % Cranck Length (mm)
Lbc = 1143; % (mm)
Lcd = 4712.396; % (mm)
Lda = 611.11; % Ground Length (mm)
theta_ab = 0:1:360; % Cranck Angle Vector (°) 
theta_ab = theta_ab*pi/180; % Cranck Angle Vector (rad)
w_ab = 10; % rad/s
alpha_ab = 0; %rad/s2


% Solución al Problema
seed = [333*pi/180,206] %Valor Semilla
for i = 1:1:361
    cranck_angle = theta_ab(i);
    [x, fval, info] = fsolve(@sol_PF,seed);
    theta_bc(i) = x(1); % Coupler Angle (rad)
    theta_cd(i) = x(2); 
    
    seed = [x(1),x(2)]; %Valor Semilla
     
    e_ab = [cos(theta_ab(i));sin(theta_ab(i))]; % Vector Columna*
    e_bc = [cos(theta_bc(i));sin(theta_bc(i))]; % Vector Columna*
    e_cd = [cos(theta_cd(i));sin(theta_cd(i))]; % Vector Columna*
    
    
    n_ab = [-sin(theta_ab(i));cos(theta_ab(i))]; % Vector Columna*
    n_bc = [-sin(theta_bc(i));cos(theta_bc(i))]; % Vector Columna*
    n_cd = [-sin(theta_cd(i));cos(theta_cd(i))]; % Vector Columna*
    
    
    % Solución del "Estado de Velocidades"
    A1 = [Lbc*n_bc, Lcd*n_cd]; % Vector Fila*
    b1 = [-w_ab*Lab*n_ab]; 
    x = A1\b1;
    w_bc(i) = x(1);
    w_cd(i) = x(2);
    
    
      % Velocidad Punto C
      Vc(:,i) = w_ab*Lab*n_ab + w_bc(i)*Lbc*n_bc;
      
      % Velocidad Punto B
      Vb(:,i) = w_cd(i)*Lcd*n_cd + w_bc(i)*Lbc*n_bc;
    
    % Solución del "Estado de Aceleraciones"
    A2 = A1; % Vector Fila*
    b2 = [-alpha_ab.*Lab.*n_ab + ((w_ab).^2).*Lab.*e_ab + ((w_bc(i)).^2).*Lbc.*e_bc + ((w_cd(i)).^2).*Lcd.*e_cd];
    x = A2\b2;
    alpha_bc(i) = x(1);
    alpha_cd(i) = x(2);
    
      % Aceleración Punto C
      Ac(:,i) = alpha_ab*Lab*n_ab - ((w_ab).^2)*Lab*e_ab + alpha_bc(i)*Lbc*n_bc - ((w_bc(i)).^2)*Lbc*e_bc;
      
      % Aceleración Punto B
      Ab(:,i) = alpha_cd(i)*Lcd*n_cd - ((w_cd(i)).^2)*Lcd*e_cd + alpha_bc(i)*Lbc*n_bc - ((w_bc(i)).^2)*Lbc*e_bc;
end

% Postprocesamiento 
figure(1)
hold on
grid on
set(gca,'FontSize',20)
xlabel('Crank Angle (°)')
ylabel('Velocidad Punto C (mm/s)')
xlim([0 360])
plot(theta_ab*180/pi,Vc(1,:),'color','r','lineWidth',4)
plot(theta_ab*180/pi,Vc(2,:),'color','b','lineWidth',4)
legend('Vx','Vy')

figure(2)
hold on
grid on
set(gca,'FontSize',20)
xlabel('Crank Angle (°)')
ylabel('Aceleración Punto C (mm/s^2)')
xlim([0 360])
plot(theta_ab*180/pi,Ac(1,:),'color','r','lineWidth',4) 
plot(theta_ab*180/pi,Ac(2,:),'color','b','lineWidth',4) 
legend('Ax','Ay')

figure(3)
hold on
grid on
set(gca,'FontSize',20)
xlabel('Crank Angle (°)')
ylabel('Velocidad Punto B (mm/s)')
xlim([0 360])
plot(theta_ab*180/pi,Vb(1,:),'color','r','lineWidth',4)
plot(theta_ab*180/pi,Vb(2,:),'color','b','lineWidth',4)
legend('Vx','Vy')

figure(4)
hold on
grid on
set(gca,'FontSize',20)
xlabel('Crank Angle (°)')
ylabel('Aceleración Punto B (mm/s^2)')
xlim([0 360])
plot(theta_ab*180/pi,Ab(1,:),'color','r','lineWidth',4) 
plot(theta_ab*180/pi,Ab(2,:),'color','b','lineWidth',4) 
legend('Ax','Ay')
