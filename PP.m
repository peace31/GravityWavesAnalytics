clc;
clear;
close all;

%% Constants for Model
L=18;
g=9.81;
H=0.75;
ro=1000;
w=1.5;
N=128;M=512;
dx=L/N;dt=0.9*dx/(g*H)^0.5;
x=0:dx:N*dx;
thita=2*pi.*x./L;
t=0:dt:(M-1)*dt;
%% Initialize of model
u_n(:,1)=5/11*(g*H)^0.5.*sin(thita./2);
h(:,1)=H*(1/3.*cos((thita-pi)./2).^16+0.2.*cos((thita-pi/2)./4).^32+(1-x./L)+691/1818);
e_ta(:,1)=h(:,1)-H;
phi(:,1)=0.01.*exp(-128.*(x./L-3/4).^2);
%% Read ASCII data
u_lin=load('u_lin.dat');
h_lin=load('h_lin.dat');
%% plot the initial value of fluid surface
figure;
subplot(3,1,1);
plot(x,u_lin(:,1));
title('u(x,0)');
subplot(3,1,2);
plot(x,h_lin(:,1));
title('h(x,0)');
subplot(3,1,3);
plot(x,phi(:,1));
title('\phi(x,0)');
%% Read the Matlab file
D=load('soln.mat');
u_nlin=D.u_nlin;
h_nlin=D.h_nlin;
phi_nlin=D.phi_nlin;
figure;
subplot(3,1,1);
plot(x,u_nlin(:,1));
title('u(x,0)');
subplot(3,1,2);
plot(x,h_nlin(:,1));
title('h(x,0)');
subplot(3,1,3);
plot(x,phi_nlin(:,1));
title('\phi(x,0)');
%% plot the evolution of the fluid depth for linear and non-lienar model
figure;
subplot(5,2,1);
plot(x,h_lin(:,1));
title('liner model');
subplot(5,2,2);
plot(x,h_nlin(:,1));
title('non-liner model');
subplot(5,2,3);
plot(x,h_lin(:,128));
subplot(5,2,4);
plot(x,h_nlin(:,128));
subplot(5,2,5);
plot(x,h_lin(:,256));
subplot(5,2,6);
plot(x,h_nlin(:,256));
subplot(5,2,7);
plot(x,h_lin(:,384));
subplot(5,2,8);
plot(x,h_nlin(:,384));
subplot(5,2,9);
plot(x,h_lin(:,512));
subplot(5,2,10);
plot(x,h_nlin(:,512));
%% visualization  by  getframe
figure;
for i =1:M
    plot(x,h_lin(:,i));
    title('liner model');
    F(i) = getframe;
end
movie(F);
figure;
for i =1:M
    plot(x,h_nlin(:,i));
    title('non-liner model');
    G(i) = getframe;
end
movie(G);
%% Creating the linear and non-linear Model
prompt = 'Enter liner/non-linear: ';
answer = input(prompt);
if(answer==1)
    for i =1:M-1
        [u2,h2,eta2,phi2] = nonlinear(u_n(:,i),h(:,i),e_ta(:,i),phi(:,i),N,dx,dt,g,H);
        u_n(:,i+1)=u2;h(:,i+1)=h2;e_ta(:,i+1)=eta2;phi(:,i+1)=phi2;
    end
    s=sum(sum(u_n-u_nlin));
    if(abs(s)<0.0001)
        disp('created model is same as saved model!');
    end
else
    for i =1:M-1
        [u2,h2,eta2] = gravity(u_n(:,i),h(:,i),e_ta(:,i),N,dx,dt,g,H);
        phi2=transport(phi(:,i),u_n(:,i),N,dx,dt);
        u_n(:,i+1)=u2;h(:,i+1)=h2;e_ta(:,i+1)=eta2;phi(:,i+1)=phi2;
    end
    s=sum(sum(u-u_lin));
    if(abs(s)<0.0001)
        disp('created model is same as saved model!');
    end
end
figure;
plot(x,phi(:,1),x,phi(:,128),x,phi(:,184));
legend('m=1','m=128','m=184');
Ke=zeros(M,1);Pe=zeros(M,1);te=zeros(M,1);mass=zeros(M,1);
for i =1:M
    %Kenetic energy
    f1=0.5.*h(:,i).*u_n(:,i).^2;
    %potential evergy
    f2=0.5*g.*e_ta(:,i).^2;
    %fluid mass 
    f3=h(:,i);
    Ke(i)=ro*w.*trapezoidal(f1,N,dx);
    Pe(i)=ro*w.*trapezoidal(f2,N,dx);
    mass(i)=ro*w.*trapezoidal(f3,N,dx);
    te(i)=Pe(i)+Ke(i);
end
figure;
subplot(1,2,1);
plot(t,mass./mass(1));
title('The Change of Fluid Mass vs time');
subplot(1,2,2);
plot(t,te./te(1));
title('The Change of Total available energy vs time');

figure;
plot(t,Ke,t,Pe,t,te);
title('The Change of Available energy vs Time');
legend('kinetic energy','potential energy','total available energy');
%% plot the evolution of the fluid depth for the created linear and non-lienar model
figure;
subplot(5,2,1);
plot(x,h_lin(:,1));
title('m=1');
subplot(5,2,2);
plot(x,h_nlin(:,1));
title('m=1');
subplot(5,2,3);
plot(x,h_lin(:,128));
title('m=128');
subplot(5,2,4);
plot(x,h_nlin(:,128));
title('m=128');
subplot(5,2,5);
plot(x,h_lin(:,256));
title('m=256');
subplot(5,2,6);
plot(x,h_nlin(:,256));
title('m=256');
subplot(5,2,7);
plot(x,h_lin(:,384));
title('m=384');
subplot(5,2,8);
plot(x,h_nlin(:,384));
title('m==384');
subplot(5,2,9);
plot(x,h_lin(:,512));
title('m=512');
subplot(5,2,10);
plot(x,h_nlin(:,512));
title('m=512');