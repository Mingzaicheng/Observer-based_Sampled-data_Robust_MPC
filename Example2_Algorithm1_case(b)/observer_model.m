function dx = observer_model(x, u, y_0, t,Lp0, hat_x_0,flag)
dx = zeros(2,1);     
A1 = [0, 1; 0.1, -0.2];
A2 = [0, 1; 0.1, 1];
AA={A1,A2};

B1=[0,1]';
B2=[0,1]';
BB=[0,1]';
CC = [0 1];
if mod(flag, 2) == 1  
    Ak = A1;  
else 
    Ak = A2;  
end
dx(1) = [1 0] * (Ak * [x(1);x(2)] + BB * u + Lp0 * (y_0 - CC * hat_x_0));
dx(2) = [0 1] * (Ak * [x(1);x(2)] + BB * u  +Lp0 * (y_0 - CC * hat_x_0));
end

