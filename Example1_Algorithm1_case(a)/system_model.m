function dx = system_model(x, u, t, dist)
dx = zeros(2,1);

beta=10;
x2bar=beta;
x2ule=-beta;

Teq=350;
CAeq=0.5;

q=100;
V=100;
ER=8750;
rho=1000;
Cp=0.239;
UA=5e4;
C_Af=1;
T_f=350;
Tceq=338;

k00=q/V*(C_Af-CAeq)/(exp(-ER/Teq)*CAeq);
DetaH=(q/V*(Teq-T_f)+UA/V/rho/Cp*(Teq-Tceq))/(k00*exp(-ER/Teq)*CAeq/rho/Cp);

fai1Bar=k00*exp(-ER/(x2bar+Teq));
fai1Ule=k00*exp(-ER/(x2ule+Teq));
fai_1_0=(fai1Bar+fai1Ule)/2;
fai2Bar=k00*(exp(-ER/(x2bar+Teq))-exp(-ER/Teq))*CAeq/x2bar;
fai2Ule=k00*(exp(-ER/(x2ule+Teq))-exp(-ER/Teq))*CAeq/x2ule;
fai_2_0=(fai2Bar+fai2Ule)/2;
g1Bar=fai1Bar-fai_1_0;
g1Ule=fai1Ule-fai_1_0;
g2Bar=fai2Bar-fai_2_0;
g2Ule=fai2Ule-fai_2_0;

g1x2=k00*exp(-ER/(x(2)+Teq))-fai_1_0;
g2x2=k00*(exp(-ER/(x(2)+Teq))-exp(-ER/Teq))*CAeq/x(2)-fai_2_0;

h1=0.5*(g1x2-g1Ule)/(g1Bar-g1Ule);
h2=0.5*(g1Bar-g1x2)/(g1Bar-g1Ule);
h3=0.5*(g2x2-g2Ule)/(g2Bar-g2Ule);
h4=0.5*(g2Bar-g2x2)/(g2Bar-g2Ule);


A1=[-q/V-fai_1_0-2*g1Bar, -fai_2_0; DetaH/rho/Cp*fai_1_0+2*DetaH/rho/Cp*g1Bar, -q/V-UA/V/rho/Cp+DetaH/rho/Cp*fai_2_0];
A2=[-q/V-fai_1_0-2*g1Ule, -fai_2_0; DetaH/rho/Cp*fai_1_0+2*DetaH/rho/Cp*g1Ule, -q/V-UA/V/rho/Cp+DetaH/rho/Cp*fai_2_0];
A3=[-q/V-fai_1_0, -fai_2_0-2*g2Bar; DetaH/rho/Cp*fai_1_0, -q/V-UA/V/rho/Cp+DetaH/rho/Cp*fai_2_0+2*DetaH/rho/Cp*g2Bar];
A4=[-q/V-fai_1_0, -fai_2_0-2*g2Ule; DetaH/rho/Cp*fai_1_0, -q/V-UA/V/rho/Cp+DetaH/rho/Cp*fai_2_0+2*DetaH/rho/Cp*g2Ule];
BB=[0, UA/V/rho/Cp]';
DD=[0.01;0.2];

Ak=h1*A1+h2*A2+h3*A3+h4*A4;


dx(1) = [1 0] * (Ak * [x(1);x(2)] + BB * u + DD * dist(t));
dx(2) = [0 1] * (Ak * [x(1);x(2)] + BB * u + DD * dist(t));
end