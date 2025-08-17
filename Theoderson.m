clc;
clear;

k=0.1:0.01:1;
m=3.3843;
Sa=0.0859;
Ia=0.0135;
rho=1.225;
c=0.2540;
b=c/2;
a=-0.5;
Kh=2818.8;
Ka=37.3;
for i=1:length(k)
    [u1(i),u2(i),g1(i),g2(i)]=U_g_method(b,k(i),Sa,Ia,m,Kh,Ka,rho,a);
end

figure(1)
hold on
plot(u1,g1,'Color','b')
plot(u2,g2,'Color','r')
yline(0,'Color','k','LineWidth',1,'Label','g=0')
xlabel('Airspeed[U]')
ylabel('g')
xlim([0,30])
box on
grid on


function C_k=Compute_Ck(k)
if k==0
    C_k=1;
    return
end
i=sqrt(-1);
H_12= besselh(1,2,k);
H_02= besselh(0,2,k);
C_k=H_12/(H_12+i*H_02);
end

function [U1,U2,g1,g2]=U_g_method(b,k,Sa,Ia,m,Kh,Ka,rho,a)
syms Omega
i=sqrt(-1);
C_k=Compute_Ck(k);
L_h=-1+2*i*C_k/k;
L_a=-0.5+i*(1+2*C_k)/k +2*C_k/k^2;
M_h=0.5;
M_a=3/8-i/k;
mu=m/(pi*rho*b^2);
x_a=Sa/(m*b);
r_a=sqrt(Ia/(m*b^2));
l_h=L_h;
omega_a2=Ka/Ia;
omega_h2=Kh/m;
l_a=L_a-(0.5+a)*L_h;
m_h=M_h+(0.5+a)*L_h;
m_a=M_a-(0.5+a)*(-L_a+M_h)-(0.5+a)^2*L_h;
K(1,1)=mu*(1-omega_h2/omega_a2*Omega)-l_h;
K(1,2)=mu*x_a-l_a;
K(2,1)=mu*x_a+m_h;
K(2,2)=mu*r_a^2*(1-Omega)+m_a;
eqn=det(K);
soln=double(solve(eqn,Omega));
X1=real(soln(1));
Y1=imag(soln(1));
X2=real(soln(2));
Y2=imag(soln(2));
U1=b*sqrt(omega_a2)/(k*sqrt(X1));
U2=b*sqrt(omega_a2)/(k*sqrt(X2));
g1=Y1/X1;
g2=Y2/X2;
end