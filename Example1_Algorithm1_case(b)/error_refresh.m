

setlmis([]);

theta_plus=lmivar(1, [1 1]);

LMI = newlmi;
lmiterm([LMI 1 1 0],0);
lmiterm([-LMI 1 1 theta_plus],1,1);

LMI = newlmi;
lmiterm([-LMI 1 1 0],1);
lmiterm([LMI 1 1 theta_plus],1,1);


for ss = 1:4
   
    LMI = newlmi;
    lmiterm([LMI 1 1 0],Pe_k0);
    lmiterm([-LMI 1 1 theta_plus],1, inv(EP_matrix{ss}));
end

dbc3=getlmis;

OP3=[0 100 1e9 10 1];

[TMIN3,VXx3_f]=feasp(dbc3,OP3);

TMIN3

if TMIN3 < 0

    CC3 = [1];
    
    [OPX3,VXx3] = mincx(dbc3, CC3, OP3, VXx3_f);
    
    theta_plus_0 = dec2mat(dbc3, VXx3, theta_plus);
    eta_plot(flag) = flag;
end


