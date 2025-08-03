
clear
Pw=100

A1 = [0, 1; 0.1, -0.2];
A2 = [0, 1; 0.1, 1];
AA={A1,A2};
B1=[0,1]';
B2=[0,1]';
BB=[0,1]';
CC = [0 1];
DD = [0.03;0.03];
alpha1=0.999;
Tsam = 0.05;
ss = 1;

Q_weight =[1,0;0,1];
R_weight = 1;
Pr = [0.5/(0.5)^2 0;0 0.5/(0.5)^2];
hat_x0 = [9;9];
e_0 = [1;1];
x_0 = hat_x0 + e_0;

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
results_H=zeros(1,20);
results_F=zeros(1,20);
Gamma0 = zeros(1,20);
eta_plot = zeros(1,20);
TMIN_results=size(1, 50);
counter=1;
T = 7;
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

    U0(flag)=F_k0*hat_x0;
    Gamma0(flag) = gamma0; 
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

        line([(iii-0.05),iii],[hat_x01(flag-1),hat_x01(flag)],[hat_x02(flag-1),hat_x02(flag)],'Color','b','linestyle',':','linewidth',1.0)
        line([(iii-0.05),iii],[x_01(flag-1),x_01(flag)],[x_02(flag-1),x_02(flag)],'Color','k','linewidth',1.0)
        hold on

    end
   
    figure (14)
    count=ones(1,201);
    kt=1:100;
    plot3(count(kt)*iii,hat_x01(flag),hat_x02(flag),'*b');
    plot3(count(kt)*iii,x_01(flag),x_02(flag),'sk'); 
    iii_values = []; 


    for ss=1:10
        
    tspan=[iii,iii+dt/2, iii+dt] 
     iii_values = [iii_values, iii];
    [~, x_r] = ode23s(@(t, x) system_model(x, U0(flag), t, dist, flag), tspan, x_0);
    x_state_plot(counter+1,:) = x_r(3,:);    
    x_0=x_r(3,:)
    
    y_0 = CC * x_0';
    [~, hat_x_r] = ode23s(@(t, x) observer_model(x, U0(flag), y_0, t, Lp0, hat_x0, flag), tspan, hat_x0);
    hat_x0_plot(counter+1,:)=hat_x_r(3,:);    
    hat_x0=hat_x_r(3,:)'
    counter=counter+1
    iii=iii+dt
    end
    
    x_0=x_0'
    
    EEr0=ellipsoid([0;0], inv(Pe_k0))
    W=ellipsoid([0], inv(Pw))
    for iii=1:2
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
legend({'Algorithm 1 for case (a),~$x_1$', 'Algorithm 1 for case (a),~$ \hat x_{1}$', 'Algorithm 1 for case (b),~$x_1$', 'Algorithm 1 for case (b),~$\hat x_{1}$'}, 'Interpreter', 'latex');
figure (11111)
plot(0:0.005:(T+0.05),x_state_plot*[0;1],'b','LineWidth',0.5)
hold on
plot(0:0.005:(T+0.05),hat_x0_plot*[0;1],'--b','LineWidth', 0.5)
hold on
xlim([0,T])
xlabel({'$t~({\rm min})$'},'Interpreter','latex');
ylabel({'$x_2,\hat x_{2}$'},'Interpreter','latex');
legend({'Algorithm 1 for case (a),~$x_2$', 'Algorithm 1 for case (a),~$ \hat x_{2}$', 'Algorithm 1 for case (b),~$x_2$', 'Algorithm 1 for case (b),~$\hat x_{2}$'}, 'Interpreter', 'latex');

figure (111)
plot(x_state_plot*[1;0],x_state_plot*[0;1],'b','LineWidth',1);
hold on
plot(hat_x0_plot*[1;0],hat_x0_plot*[0;1],'--b','LineWidth',1);
xlabel({'$x_1,\hat x_{1}$'},'Interpreter','latex');
ylabel({'$x_2,\hat x_{2}$'},'Interpreter','latex');
xlim([-0.025,0.4]);
ylim([-0.5,7.5]);
legend({'Algorithm 1 for case (a),~$x$','Algorithm 1 for case (a),~$\hat x$','Algorithm 1 for case (b),~$x$','Algorithm 1 for case (b),~$\hat x$'},'Interpreter','latex');

figure (2)
hold on
stairs(0:0.05:T, U0,'b','LineWidth', 1); 
xlabel({'$t~({\rm min})$'},'Interpreter','latex');
ylabel({'$u$'},'Interpreter','latex');
ylim([-20,5]);
legend({'Algorithm 1 for case (a)','Algorithm 1 for case (b)'},'Interpreter','latex');
box on
set(gca, 'Box', 'on') 

