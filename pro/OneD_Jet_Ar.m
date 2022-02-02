% Author:shayan dodge
% Email address:dodgeshayan@gmail.com
% 1D Jet Ar
%frequency=10KHz,v=6kv,P=1atm
clc
close all;
clear all;
format long
%% some constants
mu_0 = 1.2566370614359173e-06;
eps_0= 8.8541878176203892e-12;
C=299792458.0;
% speed of light
%% wave definition
i=1;
r=1.*10^-3;R=2.*10^-3;
amptidute=1500;
frequency=2.45.*10^9;%20.*10^9
T=1/frequency;
lambda=C*T;
omega=2*pi*frequency;
%% CPML parameters
pmlWidth=100;
pmlOrder=4;
epsR=1;
sigmaMax=1;
kappaMax=15;
alphaMax=0.24;
alphaOrder=1;
%% domain definition
domainLength=1;
totalTime=10^10*T;
numberCellsPerWavelength=20;%20
dz=lambda/numberCellsPerWavelength;
dt=dz/(2.*C);
totalTimeStep=floor(totalTime/dt);
nz=floor(domainLength/dz+2*pmlWidth);
nzp1=nz+1;
nzm1=nz-1;
zcoor=linspace (0, domainLength, (nz-2*pmlWidth)+1);

%% FDTD EM field arrarys
Ex=zeros(1,nzp1);
Hy=zeros(1,nz);
%% Coeficient for EM field updating
Cexe=1;
Cexhy=-dt/eps_0/dz;
Chyh=1;
Chyex=-dt/mu_0/dz;

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
%% Plasma parameters
pd_1z=floor(0.5./dz)+pmlWidth;
pd_2z=floor(1./dz)+pmlWidth;
zcoorp=linspace (0.5, 1, ((pd_2z-pd_1z))+1);
a=100000;c=4000000+100000;d=100000; v=([a:d:c]');
n_e=zeros(1,abs(pd_2z-pd_1z)+1);
n_e(:)=1.*10^17;
n_Ar=zeros(1,abs(pd_2z-pd_1z)+1);
n_Ar(:)=10^23;
n_Ari=zeros(1,abs(pd_2z-pd_1z)+1);
n_Ari(:)=10;%n_e(:);
k_t=1.*1.6*10^(-19);%Boltzmann constant
T_e=k_t./(1.6*10^(-19));
T_g=6.020289700000000e+02;
m_e=9.10938356*10^(-31);%electron mass (kg)
M_Ar=(39.948.*10^(-3))/(6.02214076*10^(23));%argon mass (kg)
e=1.6*10^(-19);
E_5ea=15.76*1.6*10^(-19);
E_6ea=15.76*1.6*10^(-19);
vi_5ea=(2.*v.^2+(2*E_5ea)/m_e).^(0.5);
v_2_f_0=zeros(((c-a)/d)+1,abs(pd_2z-pd_1z)+1);
v_2_f_1=zeros(((c-a)/d)+1,abs(pd_2z-pd_1z)+1);
Jx=(zeros((nz)));
for iii=1:((c-a)/d)+1
v_2_f_0(iii,:)=n_e.*(v(iii).^2.*1.*(((gamma(0.25)).^4)./(sqrt(2).*(pi.^2))).*((m_e./(12.*sqrt(2).*pi*k_t)).^(3./2)).*exp(-(((gamma(0.25)).^4)./(74.*pi.^2)).*(m_e.*v(iii).^2./(2.*k_t)).^2));
end
v_2_f_0_vi_5ea=zeros(((c-a)/d)+1,(pd_2z-pd_1z)+1);
v_2_f_0_vi_5ea(1:end-21,:)=...
      (((vi_5ea(1:end-21)-v(21:end-1))./(v(20:end-2)-v(21:end-1))).*((vi_5ea(1:end-21)-v(22:end))./(v(20:end-2)-v(22:end)))).*v_2_f_0(20:end-2,:)+...
      (((vi_5ea(1:end-21)-v(20:end-2))./(v(21:end-1)-v(20:end-2))).*((vi_5ea(1:end-21)-v(22:end))./(v(21:end-1)-v(22:end)))).*v_2_f_0(21:end-1,:)+...
      (((vi_5ea(1:end-21)-v(20:end-2))./(v(22:end)-v(20:end-2))).*((vi_5ea(1:end-21)-v(21:end-1))./(v(22:end)-v(21:end-1)))).*v_2_f_0(22:end,:);

  s_1ea=10^-20;
% load('s_1ea');
load('s_5ea');
load('s_5ea_p');
s_6ea=10^(-31).*(T_e/300).^(-4.5);

nu_i_v_5ea=n_Ar.*s_5ea.*v;
nu_i_vi_5ea=n_Ar.*s_5ea_p.*vi_5ea;
nu_Att_v_6ea=n_Ari.*s_6ea.*v;
nu_m=n_Ar.*s_1ea.*v;%6.8.*10^8;%n_Ar.*s_1ea.*v;%constant collision frequency(Hz)

k_1ea=((4.*pi./n_e).*d.*sum(s_1ea.*v(2:end).*v_2_f_0(2:end,:)));
k_5ea=((4.*pi./n_e).*d.*sum(s_5ea(2:end).*v(2:end).*v_2_f_0(2:end,:)));
k_6ea=((4.*pi./n_e).*d.*sum(s_6ea.*v(2:end).*v_2_f_0(2:end,:)));

omega_p=sqrt((max(n_e).*(e.^2))./(m_e.*eps_0))./(2.*pi)
omega
omega./omega_p

tic
%% FDTD loop
for number=1:1:floor((10.*10^-6)./dt)
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

    v_2_f_1(:,:)=(1./(2+dt.*nu_m)).*((2-dt.*nu_m).*...
    v_2_f_1(:,:)+(dt.*e./m_e).*...
    (def(v_2_f_0(:,:),v,a,c,d,nz,pd_1z,pd_2z)-2.*v_2_f_0./v).*(2.*Ex(pd_1z:pd_2z)));

    Jx(pd_1z:pd_2z)=(-4.*pi.*e./3).*d.*sum(v(1:18).*v_2_f_1(1:18,:));

    for n=1:10
dt=dt./10;
v_2_f_0(:,:)=((v_2_f_0(:,:)+...
((dt.*e)./(3.*m_e)).*((Ex(pd_1z:pd_2z)+Ex(pd_1z+1:pd_2z+1)+Ex(pd_1z-1:pd_2z-1))./3).*...
def(v_2_f_1(:,:),v,a,c,d,nz,pd_1z,pd_2z)+...
(dt).*(m_e/M_Ar).*def(nu_m.*v.^1.*v_2_f_0(:,:),v,a,c,d,nz,pd_1z,pd_2z)+...
(dt).*(T_g.*1.38.*10^-23./M_Ar).*(nu_m./v).*(def_2(v.*v_2_f_0(:,:),v,a,c,d,nz,pd_1z,pd_2z)-...
3.*def(v_2_f_0(:,:),v,a,c,d,nz,pd_1z,pd_2z))+...
4.*dt.*((vi_5ea./v)).*nu_i_vi_5ea.*v_2_f_0_vi_5ea(:,:)-...
dt.*v_2_f_0(:,:).*nu_i_v_5ea)+... 
dt.*v_2_f_0(:,:).*nu_Att_v_6ea);
dt=dt.*10;

v_2_f_0(:,:)=(n_e./(4.*pi.*d.*sum((v_2_f_0(:,:))))).*v_2_f_0(:,:);

    
nu_i_v_5ea=n_Ar.*s_5ea.*v;
nu_i_vi_5ea=n_Ar.*s_5ea_p.*vi_5ea;
nu_Att_v_6ea=n_Ari.*s_6ea.*v;
    
k_1ea=((4.*pi./n_e).*d.*sum(s_1ea.*v(2:end).*(v_2_f_0(2:end,:)+v_2_f_1(2:end,:))));
k_5ea=((4.*pi./n_e).*d.*sum(s_5ea(2:end).*v(2:end).*(v_2_f_0(2:end,:)+v_2_f_1(2:end,:))));
k_6ea=((4.*pi./n_e).*d.*sum(s_6ea.*v(2:end).*(v_2_f_0(2:end,:)+v_2_f_1(2:end,:))));
     
dn=dt*(k_5ea.*n_e.*n_Ar-k_6ea.*n_Ari.*n_e);
n_e=n_e+dn;
n_Ari=n_Ari+dt.*(k_5ea.*n_e.*n_Ar-k_6ea.*n_Ari.*n_e);
n_Ar=n_Ar-dt.*(k_5ea.*n_e.*n_Ar-k_6ea.*n_Ari.*n_e);
end

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
     Ex(pmlWidth+20)=amptidute.*sin(i*dt*omega)+Ex(pmlWidth+20);%Sine source
    %======================================================================
    % update figure
    %======================================================================

i=i+1;
gam=floor((5.*10^-7)./dt);
if mod(number,gam)==0%0.1mic
toc
plot(v_2_f_0(:,1))
A(floor(number./gam))=n_e(1);
h=figure(1)
grid on
h=plot(zcoorp,n_e)
h_2=figure(2)
grid on
h_2=plot(A)
h_3=figure(3)
grid on
h_3=plot(zcoor,Ex(1+pmlWidth:end-pmlWidth))
save(['data_2.45GHz_',num2str((number./(2.*gam))),' mics.mat']);
% print(figure(1),['n_e_2.45GHz_',num2str(floor(number./(10.*gam)))],'-dtiffn')
% savefig(figure(1),['n_e_2.45GHz_',num2str(floor(number./(10.*gam))),' mics.fig'])
% print(figure(2),['A_2.45GHz_',num2str(floor(number./(10.*gam)))],'-dtiffn')
% savefig(figure(2),['A_2.45GHz_',num2str(floor(number./(10.*gam))),' mics.fig'])
% print(figure(3),['Ex_2.45GHz_',num2str(floor(number./(10.*gam)))],'-dtiffn')
% savefig(figure(3),['Ex_2.45GHz_',num2str(floor(number./(10.*gam))),' mics.fig'])
drawnow
 end
end









