

setlmis([]);

gamma_c=lmivar(1,[1,1]);

Q_c=lmivar(1,[2 1]);
Y_c=lmivar(2,[1,2]);


Nc1=lmivar(2,[2,2]);
Nc2=lmivar(2,[2,2]);
Nc3=lmivar(2,[2,2]);
Nc4=lmivar(2,[2,2]);
Nc5=lmivar(2,[1,2]);

Ne1=lmivar(2,[2,2]);
Ne2=lmivar(2,[2,2]);
Ne3=lmivar(2,[2,2]);
Ne4=lmivar(2,[2,2]);
Ne5=lmivar(2,[1,2]);

for i=1:4
    
    LMI=newlmi;
    lmiterm([LMI 1 1 Q_c],AA{i},1,'s');
    lmiterm([LMI 1 1 Q_c],-1,1);
    lmiterm([LMI 1 1 Nc1],-1,1,'s');
    lmiterm([LMI 1 1 Q_c],AA{i},Tsam,'s');

    lmiterm([LMI 2 1 -Y_c],1,BB');
    lmiterm([LMI 2 1 Nc2],-1,1);
    lmiterm([LMI 2 1 -Nc1],1,1);
    lmiterm([LMI 2 1 Q_c],1,1);
    lmiterm([LMI 2 1 Q_c],AA{i},-Tsam);
    lmiterm([LMI 2 1 -Y_c],Tsam,BB');

    lmiterm([LMI 2 2 Nc2],1,1,'s');
    lmiterm([LMI 2 2 Q_c],-1,1);
    lmiterm([LMI 2 2 Y_c],BB,-Tsam,'s');

    lmiterm([LMI 3 1 Nc3],-1,1);
    lmiterm([LMI 3 1 -Ne1],-1,1);

    lmiterm([LMI 3 2 Nc3],1,1);
    lmiterm([LMI 3 2 -Ne2],-1,1);
    
    lmiterm([LMI 3 3 Ne3],-1,1,'s');
    lmiterm([LMI 3 3 0],AA{i}'*Pe_k0+Pe_k0*AA{i}-Pe_k0+Tsam*(Pe_k0*AA{i}+AA{i}'*Pe_k0));

    lmiterm([LMI 4 1 0],CC'*Lp0');
    lmiterm([LMI 4 1 Nc4],-1,1);
    lmiterm([LMI 4 1 -Ne1],1,1);
    lmiterm([LMI 4 1 0],Tsam*CC'*Lp0');

    lmiterm([LMI 4 2 Nc4],1,1);
    lmiterm([LMI 4 2 -Ne2],1,1);
    lmiterm([LMI 4 2 0],-Tsam*CC'*Lp0');

    lmiterm([LMI 4 3 -Ne3],1,1);
    lmiterm([LMI 4 3 Ne4],-1,1);
    lmiterm([LMI 4 3 0],Pe_k0-CC'*Lp0'*Pe_k0-Tsam*(Pe_k0*AA{i}+CC'*Lp0'*Pe_k0));

    lmiterm([LMI 4 4 Ne4],1,1,'s');
    lmiterm([LMI 4 4 0],-Pe_k0+Tsam*(Pe_k0*Lp0*CC+CC'*Lp0'*Pe_k0));

    lmiterm([LMI 5 1 Nc5],-1,1);
    
    lmiterm([LMI 5 2 Nc5],1,1);
    
    lmiterm([LMI 5 3 Ne5],-1,1);
    lmiterm([LMI 5 3 0],DD'*Pe_k0+Tsam*DD'*Pe_k0);
    
    lmiterm([LMI 5 4 Ne5],1,1);
    lmiterm([LMI 5 4 0],-Tsam*DD'*Pe_k0);
    
    lmiterm([LMI 5 5 0],-eye(1));

    lmiterm([LMI 6 1 Q_c],AA{i},Tsam);
    lmiterm([LMI 6 2 Y_c],BB,Tsam);
    lmiterm([LMI 6 4 0],Tsam*Lp0*CC);
    lmiterm([LMI 6 6 Q_c],-Tsam,1);

    lmiterm([LMI 7 3 0],Tsam*Pe_k0*AA{i});
    lmiterm([LMI 7 4 0],-Tsam*Pe_k0*Lp0*CC);
    lmiterm([LMI 7 5 0],Tsam*Pe_k0*DD);
    lmiterm([LMI 7 7 0],-Tsam*Pe_k0);

    lmiterm([LMI 8 1 Q_c],Q_weight^0.5,1);
    lmiterm([LMI 8 8 gamma_c],-1,eye(1));

    lmiterm([LMI 9 2 Y_c],R_weight^0.5,1);
    lmiterm([LMI 9 9 gamma_c],-1,eye(1));
    
end

for i=1:4
    
    LMI=newlmi;
    lmiterm([LMI 1 1 Q_c],AA{i},1,'s');
    lmiterm([LMI 1 1 Q_c],-1,1);
    lmiterm([LMI 1 1 Nc1],-1,1,'s');

    lmiterm([LMI 2 1 -Y_c],1,BB');
    lmiterm([LMI 2 1 Nc2],-1,1);
    lmiterm([LMI 2 1 -Nc1],1,1);
    lmiterm([LMI 2 1 Q_c],1,1);

    lmiterm([LMI 2 2 Nc2],1,1,'s');
    lmiterm([LMI 2 2 Q_c],-1,1);

    lmiterm([LMI 3 1 Nc3],-1,1);
    lmiterm([LMI 3 1 -Ne1],-1,1);

    lmiterm([LMI 3 2 Nc3],1,1);
    lmiterm([LMI 3 2 -Ne2],-1,1);

    lmiterm([LMI 3 3 Ne3],-1,1,'s');
    lmiterm([LMI 3 3 0],AA{i}'*Pe_k0+Pe_k0*AA{i}-Pe_k0);

    lmiterm([LMI 4 1 0],CC'*Lp0');
    lmiterm([LMI 4 1 Nc4],-1,1);
    lmiterm([LMI 4 1 -Ne1],1,1);

    lmiterm([LMI 4 2 Nc4],1,1);
    lmiterm([LMI 4 2 -Ne2],1,1);

    lmiterm([LMI 4 3 -Ne3],1,1);
    lmiterm([LMI 4 3 Ne4],-1,1);
    lmiterm([LMI 4 3 0],Pe_k0-CC'*Lp0'*Pe_k0);

    lmiterm([LMI 4 4 Ne4],1,1,'s');
    lmiterm([LMI 4 4 0],-Pe_k0);

    lmiterm([LMI 5 1 Nc5],-1,1);
    
    lmiterm([LMI 5 2 Nc5],1,1);
    
    lmiterm([LMI 5 3 Ne5],-1,1);
    lmiterm([LMI 5 3 0],DD'*Pe_k0);
    
    lmiterm([LMI 5 4 Ne5],1,1);
    
    lmiterm([LMI 5 5 0],-eye(1));

    lmiterm([LMI 6 1 Q_c],AA{i},Tsam);
    lmiterm([LMI 6 2 Y_c],BB,Tsam);
    lmiterm([LMI 6 4 0],Tsam*Lp0*CC);
    lmiterm([LMI 6 6 Q_c],-Tsam,1);

    lmiterm([LMI 7 3 0],Tsam*Pe_k0*AA{i});
    lmiterm([LMI 7 4 0],-Tsam*Pe_k0*Lp0*CC);
    lmiterm([LMI 7 5 0],Tsam*Pe_k0*DD);
    lmiterm([LMI 7 7 0],-Tsam*Pe_k0);

    lmiterm([LMI 8 1 Q_c],Q_weight^0.5,1);
    lmiterm([LMI 8 8 gamma_c],-1,eye(1));

    lmiterm([LMI 9 2 Y_c],R_weight^0.5,1);
    lmiterm([LMI 9 9 gamma_c],-1,eye(1));

    lmiterm([LMI 10 1 -Nc1],Tsam,1);
    lmiterm([LMI 10 2 -Nc2],Tsam,1);
    lmiterm([LMI 10 3 -Nc3],Tsam,1);
    lmiterm([LMI 10 4 -Nc4],Tsam,1);
    lmiterm([LMI 10 5 -Nc5],Tsam,1);
    lmiterm([LMI 10 10 Q_c],-Tsam,1);

    lmiterm([LMI 11 1 -Ne1],1,Tsam);
    lmiterm([LMI 11 2 -Ne2],1,Tsam);
    lmiterm([LMI 11 3 -Ne3],1,Tsam);
    lmiterm([LMI 11 4 -Ne4],1,Tsam);
    lmiterm([LMI 11 5 -Ne5],1,Tsam);
    lmiterm([LMI 11 11 0],-Tsam*Pe_k0);
    
end


LMI=newlmi;
lmiterm([-LMI 1 1 0],1);
lmiterm([-LMI 2 1 0],hat_x0);
lmiterm([-LMI 2 2 Q_c],1,1);

LMI=newlmi;
lmiterm([-LMI 1 1 Q_c],1,1);
lmiterm([-LMI 2 1 Y_c],1,1);
lmiterm([-LMI 2 2 0],100);


dbc2=getlmis;
OP2=[0 100 1e9 10 1];
[TMIN2,VXx2_f]=feasp(dbc2,OP2);
