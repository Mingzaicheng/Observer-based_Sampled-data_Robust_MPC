
clc
clear

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

ACont1=[-q/V-fai_1_0-2*g1Bar, -fai_2_0; DetaH/rho/Cp*fai_1_0+2*DetaH/rho/Cp*g1Bar, -q/V-UA/V/rho/Cp+DetaH/rho/Cp*fai_2_0];
ACont2=[-q/V-fai_1_0-2*g1Ule, -fai_2_0; DetaH/rho/Cp*fai_1_0+2*DetaH/rho/Cp*g1Ule, -q/V-UA/V/rho/Cp+DetaH/rho/Cp*fai_2_0];
ACont3=[-q/V-fai_1_0, -fai_2_0-2*g2Bar; DetaH/rho/Cp*fai_1_0, -q/V-UA/V/rho/Cp+DetaH/rho/Cp*fai_2_0+2*DetaH/rho/Cp*g2Bar];
ACont4=[-q/V-fai_1_0, -fai_2_0-2*g2Ule; DetaH/rho/Cp*fai_1_0, -q/V-UA/V/rho/Cp+DetaH/rho/Cp*fai_2_0+2*DetaH/rho/Cp*g2Ule];
BCont=[0, UA/V/rho/Cp]';

Tsam = 0.05;
ss = 1;
AA = {ss*ACont1,ss*ACont2,ss*ACont3,ss*ACont4};
BB = BCont;
CC = [0 1];
DD = [0.01;0.2];
Pw=1;

Q_weight =[25,0;0,25];
R_weight = 1;
Pr = [0.5/(0.12)^2 0;0 0.5/(1.2)^2];
hat_x0 = [0.32;6.4]; 
 e_0 = [0.12;1.2];
x_0 = hat_x0 + [0.12;1.2];

Input = 1^2;
dist = @(t) sin(t);   

observer_off
TMIN1
if TMIN1<0
    CC1(1)=1;
    n_sub1=decnbr(dbc1);
    CC1(n_sub1)=0;
    [OPX1,VXx1]=mincx(dbc1,CC1,OP1,VXx1_f);
    Ye_k0=dec2mat(dbc1,VXx1,Ye_k)
    Pe_k0=dec2mat(dbc1,VXx1,Pe_k)
    Lp0 = inv(Pe_k0)*Ye_k0
end
 disp(Lp0)
Elliposid_Me0=ellipsoid([0;0], inv(Pr));
figure (10);
plot(Elliposid_Me0,'r')     
hold on
Elliposid_Pe0=ellipsoid([0;0], inv(Pe_k0));
figure (10);
plot(Elliposid_Pe0,'k')    
hold on


x_state_plot = zeros(211,2);
x_state_plot(1,:) = x_0';
hat_x0_plot = zeros(211,2);
hat_x0_plot(1,:) = hat_x0;
hat_x01 = zeros(1,10);
hat_x02 = zeros(1,10);
x_01 = zeros(1,10);
x_02 = zeros(1,10);
U0 = zeros(1,20);
eta_plot = zeros(1,20);
TMIN_results=size(1, 50);
counter=1;
T = 2;
dt = 0.005;
flag = 1;

for iii=0:0.05:T
    
    optimization_on
    
    TMIN2
    
    if TMIN2<0
        TMIN_results(flag) = TMIN2;
        CC2(1)=1;
        n_sub2=decnbr(dbc2);
        CC2(n_sub2)=0;
        [OPX2,VXx2]=mincx(dbc2,CC2,OP2,VXx2_f);
        gamma0=dec2mat(dbc2,VXx2,gamma_c);
        Q_c0=dec2mat(dbc2,VXx2,Q_c);
        Y_c0=dec2mat(dbc2,VXx2,Y_c);
        F_k0=Y_c0*inv(Q_c0);
    end
   
     disp(F_k0);
     disp(Y_c0);
     results_F(flag)=F_k0*hat_x0;
    U0(flag)=sign(F_k0*hat_x0)*min(abs(F_k0*hat_x0),10);
    Gamma0(flag) = gamma0; 
    figure (10)
    hat_x01(flag) = [1 0]*hat_x0;
    hat_x02(flag) = [0 1]*hat_x0;
    x_01(flag) = [1 0]*x_0;
    x_02(flag) = [0 1]*x_0;
    
    Elliposid_P1=ellipsoid(hat_x0, inv(Pe_k0));
    figure (10);
    plot(Elliposid_P1,'k');
    
    figure_info=findall(gcf,'type','line');
    xc=get(figure_info,'xdata');
    yc=get(figure_info,'ydata');
    XX=xc(2*1,1);
    YY=yc(2*1,1);
    
    Data=[XX{1,:};YY{1,:}];
    
    [z, a, b, alpha_c] = fitellipse(Data, 'linear', 'constraint', 'trace');
    
    figure (11)
    plotellipse(z, a, b, alpha_c,'k')
    
    figure_info=findall(gcf,'type','line');
    xxc=get(figure_info,'xdata');
    yyc=get(figure_info,'ydata');
    DDXX=xxc;
    DDYY=yyc;
    
    figure (14)
    count=ones(1,201);
    kt=1:100;
    plot3(count(kt)*iii,DDXX(kt),DDYY(kt),'r');
    hold on
    
    if flag>1 

        line([(iii-0.05),iii],[hat_x01(flag-1),hat_x01(flag)],[hat_x02(flag-1),hat_x02(flag)],'Color','b','linestyle',':','linewidth',1.0)%蓝色虚线
        line([(iii-0.05),iii],[x_01(flag-1),x_01(flag)],[x_02(flag-1),x_02(flag)],'Color','k','linewidth',1.0) %黑色实现
        hold on

    end

    figure (14)
    count=ones(1,201);
    kt=1:100;
    plot3(count(kt)*iii,hat_x01(flag),hat_x02(flag),'*b'); 
    plot3(count(kt)*iii,x_01(flag),x_02(flag),'sk'); 
    
    
    for ss=1:10
    tspan=[iii,iii+dt/2, iii+dt]    
    [~, x_r] = ode23s(@(t, x) system_model(x, U0(flag), t, dist), tspan, x_0);
    x_state_plot(counter+1,:) = x_r(3,:);    
    x_0=x_r(3,:)
    
    y_0 = CC * x_0';
    [~, hat_x_r] = ode23s(@(t, x) observer_model(x, U0(flag),y_0, Lp0, hat_x0), tspan, hat_x0);
    hat_x0_plot(counter+1,:)=hat_x_r(3,:);    
    hat_x0=hat_x_r(3,:)'
    counter=counter+1
    iii=iii+dt
    end
    
    x_0=x_0'
    

    EEr0=ellipsoid([0;0], inv(Pe_k0))
    W=ellipsoid([0], inv(Pw))
    for iii=1:4
        lsys1 = linsys(AA{iii}-Lp0*CC, DD, W);
        L0=[1 0; 0 1];
        TS=[(flag - 1)*Tsam, flag*Tsam] ;
        rs1 = reach(lsys1, EEr0, L0, TS);  
        E_Tsam=cut(rs1, flag*Tsam)
        
        ee=get_ea(E_Tsam)
        [Ce_1,Pe_1]=double(ee(1,1))
        [Ce_2,Pe_2]=double(ee(2,1))
        
        EPe_1=ellipsoid([0;0], Pe_1)
        EPe_2=ellipsoid([0;0], Pe_2)
        EP{iii}=[EPe_1,EPe_2]
        
        EP_matrix{iii}=double(ellunion_ea(EP{iii}));
        EPe{iii}=ellipsoid([0;0], EP_matrix{iii})
        
        hold on
    end

    error_refresh
    
    Pe_k0=1/theta_plus_0*Pe_k0;
    flag = flag + 1
    
end

figure (1)
plot(0:0.005:(T+0.05), x_state_plot*[1;0],'b','LineWidth', 0.5); 
hold on
plot(0:0.005:(T+0.05), hat_x0_plot*[1;0],'--b','LineWidth', 0.5); 
xlim([0,T])
xlabel({'$t~({\rm min})$'},'Interpreter','latex');
ylabel({'$x_1,\hat x_{1}$'},'Interpreter','latex');
legend({'Algorithm 1 for case (a),~$x_1$', 'Algorithm 1 for case (a),~$ \hat x_{1}$'}, 'Interpreter', 'latex');
figure (11111)
plot(0:0.005:(T+0.05),x_state_plot*[0;1],'b','LineWidth',0.5)
hold on
plot(0:0.005:(T+0.05),hat_x0_plot*[0;1],'--b','LineWidth', 0.5)
hold on
xlim([0,T])
xlabel({'$t~({\rm min})$'},'Interpreter','latex');
ylabel({'$x_2,\hat x_{2}$'},'Interpreter','latex');
legend({'Algorithm 1 for case (a),~$x_2$', 'Algorithm 1 for case (a),~$ \hat x_{2}$'}, 'Interpreter', 'latex');

figure (111)
plot(x_state_plot*[1;0],x_state_plot*[0;1],'b','LineWidth',1);
hold on
plot(hat_x0_plot*[1;0],hat_x0_plot*[0;1],'--b','LineWidth',1);
xlabel({'$x_1,\hat x_{1}$'},'Interpreter','latex');
ylabel({'$x_2,\hat x_{2}$'},'Interpreter','latex');
xlim([-0.025,0.4]);
ylim([-0.5,7.5]);
legend({'Algorithm 1 for case (a),~$x$','Algorithm 1 for case (a),~$\hat x$'},'Interpreter','latex');

figure (2)
hold on
stairs(0:0.05:T, U0,'--b','LineWidth', 1); 
xlabel({'$t~({\rm min})$'},'Interpreter','latex');
ylabel({'$u$'},'Interpreter','latex');
ylim([-10,1]);
legend({'Algorithm 1 for case (a)'},'Interpreter','latex');
box on
set(gca, 'Box', 'on') 

figure (14)
grid on;
set(gca,'GridLineStyle',':','GridColor','r','GridAlpha',0.3);
xlabel({'$t~({\rm min})$'},'Interpreter','latex');
xlim([0,T])
ylabel({'$x_1,\hat x_{1}$'},'Interpreter','latex');
ylim([-0.2,0.5])
zlabel({'$x_2,\hat x_{2}$'},'Interpreter','latex');
zlim([-2,12]);