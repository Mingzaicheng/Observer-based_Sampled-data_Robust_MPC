setlmis([]);

gamma_e=lmivar(1,[1 1]);
Pe_k=lmivar(1,[2 1]);

Ye_k=lmivar(2,[2 1]);

N1=lmivar(2,[2,2]);
N2=lmivar(2,[2,2]);
N3=lmivar(2,[1,2]);


for i=1:2
    LMI=newlmi;
    lmiterm([LMI 1 1 Pe_k],1,AA{i},'s');
    lmiterm([LMI 1 1 N1],-1,1,'s');
    lmiterm([LMI 1 1 Pe_k],-1,1);
    lmiterm([LMI 1 1 Pe_k],Tsam,AA{i},'s');

    lmiterm([LMI 2 1 -Ye_k],CC',-1);
    lmiterm([LMI 2 1 Pe_k],1,1);
    lmiterm([LMI 2 1 N2],-1,1);
    lmiterm([LMI 2 1 -N1],1,1);
    lmiterm([LMI 2 1 Pe_k],-Tsam,AA{i});
    lmiterm([LMI 2 1 -Ye_k],-Tsam*CC',1);

    lmiterm([LMI 2 2 N2],1,1,'s');
    lmiterm([LMI 2 2 Pe_k],-1,1);
    lmiterm([LMI 2 2 Ye_k],Tsam,CC,'s');
    
    lmiterm([LMI 3 1 Pe_k],DD',1);
    lmiterm([LMI 3 1 N3],-1,1);
    lmiterm([LMI 3 1 Pe_k],DD',Tsam);
    
    lmiterm([LMI 3 2 N3],1,1);
    lmiterm([LMI 3 2 Pe_k],DD',-Tsam);
    
    lmiterm([LMI 3 3 0],-Pw);

    lmiterm([LMI 4 1 Pe_k],Tsam, AA{i});
    lmiterm([LMI 4 2 Ye_k],-Tsam,CC);
    lmiterm([LMI 4 3 Pe_k],Tsam,DD);
    lmiterm([LMI 4 4 Pe_k],Tsam,-1);
    
    lmiterm([LMI 5 1 0],eye(2));
    lmiterm([LMI 5 5 gamma_e],1,-eye(2));
end

for i=1:2
    LMI=newlmi;
    lmiterm([LMI 1 1 Pe_k],1,AA{i},'s');
    lmiterm([LMI 1 1 N1],-1,1,'s');
    lmiterm([LMI 1 1 Pe_k],-1,1);

    lmiterm([LMI 2 1 -Ye_k],CC',-1);
    lmiterm([LMI 2 1 Pe_k],1,1);
    lmiterm([LMI 2 1 N2],-1,1);
    lmiterm([LMI 2 1 -N1],1,1);

    lmiterm([LMI 2 2 N2],1,1,'s');
    lmiterm([LMI 2 2 Pe_k],-1,1);
    
    lmiterm([LMI 3 1 Pe_k],DD',1);
    lmiterm([LMI 3 1 N3],-1,1);
    
    lmiterm([LMI 3 2 N3],1,1);
    
    lmiterm([LMI 3 3 0],-Pw);

    lmiterm([LMI 4 1 Pe_k],Tsam, AA{i});
    lmiterm([LMI 4 2 Ye_k],-Tsam,CC);
    lmiterm([LMI 4 3 Pe_k],Tsam,DD);
    lmiterm([LMI 4 4 Pe_k],Tsam,-1);
    
    lmiterm([LMI 5 1 -N1],Tsam,1);
    lmiterm([LMI 5 2 -N2],Tsam,1);
    lmiterm([LMI 5 3 -N3],Tsam,1);
    lmiterm([LMI 5 5 Pe_k],Tsam,-1);  
    
    lmiterm([LMI 6 1 0],eye(2));
    lmiterm([LMI 6 6 gamma_e],1,-eye(2));

end

LMI=newlmi;
lmiterm([-LMI 1 1 0],Pr);
lmiterm([LMI 1 1 Pe_k],1,1);

dbc1=getlmis;
OP1=[0 100 1e9 10 1];
[TMIN1,VXx1_f]=feasp(dbc1,OP1);
