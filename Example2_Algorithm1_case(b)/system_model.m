function dx = system_model(x, u, t, dist,flag)
dx = zeros(2,1);
A1 = [0, 1; 0.1, -0.2];
A2 = [0, 1; 0.1, 1];
AA={A1,A2};
B1=[0,1]';
B2=[0,1]';
BB=[0,1]';
CC = [0 1];
DD = [0.03;0.03];
if mod(flag, 2) == 1  
    Ak = A1;  
else  
    Ak = A2;  
end

dx(1) = [1 0] * (Ak * [x(1);x(2)] + BB * u + DD * dist(t));
dx(2) = [0 1] * (Ak * [x(1);x(2)] + BB * u + DD * dist(t));

end