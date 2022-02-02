
% FDTD 1D in plasma with CPML
% Author:shayan dodge
% Email address:dodgeshayan@gmail.com
%% 
clc
close all;
clear all;
%% some constants
mu_0 = 1.2566370614359173e-06;   %
eps_0= 8.8541878176203892e-12;   %
C=299792458.0;      % speed of light
%% wave definition
amptidute=20;
frequency=5*10^(10);
T=1/frequency;
lambda=C*T;
t0=4*T;
omega=2*pi*frequency;
%% CPML parameters
pmlWidth=40;
pmlOrder=4;
epsR=1;
sigmaMax=1;
kappaMax=15;
alphaMax=0.24;
alphaOrder=1;
%% domain definition
domainLength=30*lambda;
totalTime=1000*T;
numberCellsPerWavelength=20;
dz=lambda/numberCellsPerWavelength;
dt=dz/2/C;
totalTimeStep=floor(totalTime/dt);
nz=floor(domainLength/dz+2*pmlWidth);
nzp1=nz+1;
nzm1=nz-1;
ksource = pmlWidth+0;
%% FDTD EM field arrarys
Ex=zeros(1,nzp1);
Hy=zeros(1,nz);
%% Coeficient for EM field updating
Cexe=1;
Cexhy=-dt/eps_0/dz;
Chyh=1;
Chyex=-dt/mu_0/dz;
%% Plasma parameters and functions
N=1.86*10^18;
k_t=1.6*10^(-19);%Boltzmann constant
nu=6.8*10^8;%constant collision frequency(Hz)
m=9.10938356*10^(-31);%electron mass (kg)
e=1.6*10^(-19);
syms v
f_0=N.*((m/(2*pi*k_t))^(1.5))*exp((-m*(v^2))/(2*k_t));
R_f_0=diff(f_0);
f_1=sym(zeros(1,(nz)+1));
F_1=zeros(1,(nz)+1);
Jx=zeros(1,(nz)+1);
%% CPML arrays
Psi_Ezy_zn=zeros(1,pmlWidth);
Psi_Ezy_zp=zeros(1,pmlWidth);
Psi_Hzx_zn=zeros(1,pmlWidth);
Psi_Hzx_zp=zeros(1,pmlWidth);
%% initial PML update coeficients
sigmaOpt=sigmaMax*(pmlOrder+1)/(sqrt(epsR)*150*pi*dz);
% zn side
rho_e=((pmlWidth:-1:1)-0.75)/pmlWidth;
rho_m=((pmlWidth:-1:1)-0.25)/pmlWidth;
sigma_e=sigmaOpt*abs(rho_e).^pmlOrder;
sigma_m=sigmaOpt*abs(rho_m).^pmlOrder;
kappa_e=1+(kappaMax-1)*abs(rho_e).^pmlOrder;
kappa_m=1+(kappaMax-1)*abs(rho_m).^pmlOrder;
alpha_e=alphaMax*abs(rho_e).^pmlOrder;
alpha_m=alphaMax*abs(rho_m).^pmlOrder;
cpml_b_e_n=exp((-dt/eps_0)...
    *((sigma_e./kappa_e)+alpha_e));
cpml_a_e_n=1/dz*(cpml_b_e_n-1).*sigma_e ...
    ./(kappa_e.*(sigma_e+kappa_e.*alpha_e));
cpml_b_m_n=exp((-dt/eps_0)...
    *((sigma_m./kappa_m)+alpha_m));
cpml_a_m_n=1/dz*(cpml_b_m_n-1).*sigma_m ...
    ./(kappa_m.*(sigma_m+kappa_m.*alpha_m));
b=0;%sigma_e*dt/2/eps_0;
Cexe_zn=(1-b)./(1+b);
CPsi_Ezy_zn=-dt./(1+b)/eps_0;
Cexhy_zn=CPsi_Ezy_zn./kappa_e/dz;
b=0;%sigma_m*dt/2/mu_0;
Chyh_zn=(1-b)./(1+b);
CPsi_Hzx_zn=-dt./(1+b)/mu_0;
Chyex_zn=CPsi_Hzx_zn./kappa_m/dz;
% zp side
rho_e=((1:1:pmlWidth)-0.75)/pmlWidth;
rho_m=((1:1:pmlWidth)-0.25)/pmlWidth;
sigma_e=sigmaOpt*abs(rho_e).^pmlOrder;
sigma_m=sigmaOpt*abs(rho_m).^pmlOrder;
kappa_e=1+(kappaMax-1)*abs(rho_e).^pmlOrder;
kappa_m=1+(kappaMax-1)*abs(rho_m).^pmlOrder;
alpha_e=alphaMax*abs(rho_e).^pmlOrder;
alpha_m=alphaMax*abs(rho_m).^pmlOrder;
cpml_b_e_p=exp((-dt/eps_0)...
    *((sigma_e./kappa_e)+alpha_e));
cpml_a_e_p=1/dz*(cpml_b_e_p-1).*sigma_e ...
    ./(kappa_e.*(sigma_e+kappa_e.*alpha_e));
cpml_b_m_p=exp((-dt/eps_0)...
    *((sigma_m./kappa_m)+alpha_m));
cpml_a_m_p=1/dz*(cpml_b_m_p-1).*sigma_m ...
    ./(kappa_m.*(sigma_m+kappa_m.*alpha_m));
b=0;%sigma_e*dt/2/eps_0;
Cexe_zp=(1-b)./(1+b);
CPsi_Ezy_zp=-dt./(1+b)/eps_0;
Cexhy_zp=CPsi_Ezy_zp./kappa_e/dz;
b=0;%sigma_m*dt/2/mu_0;
Chyh_zp=(1-b)./(1+b);
CPsi_Hzx_zp=-dt./(1+b)/mu_0;
Chyex_zp=CPsi_Hzx_zp./kappa_m/dz;
%% Initial Plot
h_1=figure;
h_2=plot(Ex);
% set(gca,'ylim',[-1 1]*amptidute);
set(gca,'xlim',[pmlWidth (nz-pmlWidth)]);
grid on;
title(gca,strcat('E_x'));
h_3=figure;
h_4=plot(Hy);
% set(gca,'ylim',[-1 1]);
set(gca,'xlim',[pmlWidth (nz-pmlWidth)]);
grid on;
title(gca,strcat('H_y'));
h_5=figure;
h_6=plot(F_1);
% set(gca,'ylim',[-1.5 1.5]);
set(gca,'xlim',[pmlWidth (nz-pmlWidth)]);
grid on;
title(gca,strcat('F_1'));
h_7=figure;
h_8=plot(Jx);
% set(gca,'ylim',[-1.5 1.5]);
set(gca,'xlim',[pmlWidth (nz-pmlWidth)]);
grid on;
title(gca,strcat('Jx'));
%% FDTD loop
for n=1:totalTimeStep
    %======================================================================
    %update Hy
    %======================================================================
    Psi_Hzx_zn=cpml_b_m_n.*Psi_Hzx_zn+cpml_a_m_n.*(Ex(2:pmlWidth+1)-...
        Ex(1:pmlWidth));
    Psi_Hzx_zp=cpml_b_m_p.*Psi_Hzx_zp+cpml_a_m_p.*(Ex(nzp1-pmlWidth+...
        1:nzp1)-Ex(nzp1-pmlWidth:nz));
    Hy(1:pmlWidth)=Chyh_zn.*Hy(1:pmlWidth)+Chyex_zn.*(Ex(2:pmlWidth+1)-...
        Ex(1:pmlWidth))+CPsi_Hzx_zn.*Psi_Hzx_zn;
    Hy(nz-pmlWidth+1:nz)=Chyh_zp.*Hy(nz-pmlWidth+1:nz)+Chyex_zp.*...
        (Ex(nzp1-pmlWidth+1:nzp1)-Ex(nzp1-pmlWidth:nz))+CPsi_Hzx_zp.*...
        Psi_Hzx_zp;
    % non pml region
    Hy(pmlWidth+1:nz-pmlWidth)=Chyh* Hy(pmlWidth+1:nz-pmlWidth)+Chyex*...
        ( Ex(pmlWidth+2:nzp1-pmlWidth)-Ex(pmlWidth+1:nz-pmlWidth)); 
    f_1(351:end)=(1./(2+dt.*nu)).*((2-dt.*nu).*f_1(351:end)+((dt*e)/m).*...
        R_f_0.*(Ex(351:nz+1)+Ex(350:nz)));
    Jx(351:end)=(-4*pi*e/3)*int((v^3)*f_1(351:end),0,inf);
    %===========================
    % update Ex
    %===========================
    Psi_Ezy_zn=cpml_b_e_n.*Psi_Ezy_zn+cpml_a_e_n.*(Hy(2:pmlWidth+1)-...
        Hy(1:pmlWidth));
    Psi_Ezy_zp=cpml_b_e_p.*Psi_Ezy_zp+cpml_a_e_p.*(Hy(nzp1-pmlWidth:nz)-...
        Hy(nz-pmlWidth:nzm1));
    Ex(2:pmlWidth+1)=Cexe_zn.*Ex(2:pmlWidth+1)+Cexhy_zn.*(Hy(2:pmlWidth+1)...
        -Hy(1:pmlWidth))+CPsi_Ezy_zn.*Psi_Ezy_zn-(dt/(2*eps_0)).*...
        (Jx(1:pmlWidth)+Jx(2:pmlWidth+1));
    Ex(nzp1-pmlWidth:nz)=Cexe_zp.*Ex(nzp1-pmlWidth:nz)+Cexhy_zp.*...
        (Hy(nzp1-pmlWidth:nz)-Hy(nz-pmlWidth:nzm1))+CPsi_Ezy_zp.*...
        Psi_Ezy_zp-(dt/(2*eps_0)).*(Jx((nz)-pmlWidth+1:(nz))+Jx((nz)...
        -pmlWidth:(nz)-1));
    % non pml region
    Ex(pmlWidth+2:nz-pmlWidth)=Cexe* Ex(pmlWidth+2:nz-pmlWidth)+Cexhy*...
        ( Hy(pmlWidth+2:nz-pmlWidth)-Hy(pmlWidth+1:nzm1-pmlWidth))-...
        (dt/(2*eps_0)).*(Jx(pmlWidth+2:nz-pmlWidth)+Jx(pmlWidth+1:nz-1-...
        pmlWidth)); 
    
    %======================================================================
    % update source
    %======================================================================
     Ex(ksource)=Ex(ksource)-amptidute *sin((n * dt - t0) * omega);%Sine source
    %======================================================================
    % update figure
    %======================================================================
    if mod(n,100)==0||n==totalTimeStep
        set(h_2,'YData',Ex);
        pause(0.2);
        set(h_4,'YData',Hy);
        pause(0.2);
        F_1=subs(f_1,v,10^3);
        set(h_6,'YData',F_1);
        pause(0.2);  
        set(h_8,'YData',Jx);     
        pause(0.2);
        
        saveas(h_1,'Eplasmaprime_50 ','fig') 
        saveas(h_3,'Hplasmaprime_50 ','fig')
        saveas(h_5,'F_1plasmaprime_50 ','fig') 
        saveas(h_7,'Jxplasmaprime_50 ','fig')        
        save('plasma_prime_all_variable_50')
    end
    disp(n)
end
