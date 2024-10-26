clear;clc;
load("grid_a.mat");load("grid_v.mat");load("grid_x.mat");

% a = a(2:end-1); v = v(2:end-1); x = x(2:end-1); %去掉两端点，如果不想去，请注释
y = v.^2;
n = numel(x) - 1;
h =diff(x)';
delta =diff(y)';
dy0 = 2*a(1);
dy1 = 2*a(end);
ddy0 = 0;
ddy1 = 0;

%% 三次样条(两端四次)
alpha = zeros(n-1,1);
alpha(1) = -3*h(1)^3;
for i = 2:n-1
    alpha(i) = h(i);
end
alpha_m = [diag(alpha),zeros(n-1,2)];

beta = zeros(n-1,1);
for i = 1:n-2
    beta(i) = h(i+1);
end
beta(n-1) = -3*h(n)^3;
beta_m = [zeros(n-1,2),diag(beta)];

gamma = zeros(n-1,1);
for i = 1:n-1
    gamma(i) = 2*(h(i)+h(i+1));
end
gamma_m = [zeros(n-1,1),diag(gamma),zeros(n-1,1)];
A = alpha_m + beta_m + gamma_m;
A = [[-3*h(1)^3 , h(1), zeros(1,n-1)];A;[zeros(1,n-1),h(n),-3*h(n)^3]];

theta = zeros(n-1,1);
theta(1) = 3*(delta(2)/h(2) - delta(1)/h(1)) - h(1)*ddy0/2;
for i = 2:n-2
    theta(i) = 3*(delta(i+1)/h(i+1) - delta(i)/h(i));
end
theta(n-1) = 3*(delta(n)/h(n) - delta(n-1)/h(n-1)) - h(n)*ddy1/2;
theta = [3*(delta(1)/h(1) - dy0) - ddy0*h(1); theta ; 3*(dy1 - delta(n)/h(n)) - ddy1*h(n)];
M = A\theta;


e1 = M(1);
en = M(n+1);

c(1) = ddy0/2;
c(2:n) = M(2:n);

d(1) = 1/3*(c(2) - c(1))/h(1) - 2*e1*h(1);
d(2:n-1) = 1/3*(c(3:n) - c(2:n-1))./h(2:n-1);
d(n) = 1/3*(ddy1/2 - c(n))/h(n) - 2*en*h(n);%二阶连续条件

b(1) = delta(1)/h(1) - 1/3*(2*c(1)+c(2))*h(1) + e1*h(1)^3;
b(2:n-1) = delta(2:n-1)./h(2:n-1) - 1/3*(2*c(2:n-1)+c(3:n)).*h(2:n-1);
b(n) = delta(n)./h(n) - 1/3*(2*c(n)+ddy1/2).*h(n) + en*h(n)^3;%零阶重合条件

N = 20;
for k = 1:n
    ds = h(k)/N;
    if k == 1 
        for i = 1:N
            s = x(k) + (i-1)*ds;
            profile_position((k-1)*N + i) = s;
            profile_R((k-1)*N + i) = y(k) + b(k)*(s-x(k)) + c(k)*(s-x(k))^2 ...
                                     + d(k)*(s-x(k))^3 + e1*(s-x(k))^4;
            profile_Rd1((k-1)*N + i) = b(k) + 2*c(k)*(s-x(k)) ...
                                 + 3*d(k)*(s-x(k))^2 + 4*e1*(s-x(k))^3;
            profile_Rd2((k-1)*N + i) = 2*c(k) ...
                                 + 6*d(k)*(s-x(k)) + 12*e1*(s-x(k))^2;
        end
        continue;
    end
     if k == n 
        for i = 1:N
            s = x(k) + (i-1)*ds;
            profile_position((k-1)*N + i) = s;
            profile_R((k-1)*N + i) = y(k) + b(k)*(s-x(k)) + c(k)*(s-x(k))^2 ...
                                     + d(k)*(s-x(k))^3 + en*(s-x(k))^4;
            profile_Rd1((k-1)*N + i) = b(k) + 2*c(k)*(s-x(k)) ...
                                 + 3*d(k)*(s-x(k))^2 + 4*en*(s-x(k))^3;
            profile_Rd2((k-1)*N + i) = 2*c(k) ...
                                 + 6*d(k)*(s-x(k)) + 12*en*(s-x(k))^2;
        end
        break;
    end
    for i = 1:N
        s = x(k) + (i-1)*ds;
        profile_position((k-1)*N + i) = s;
        profile_R((k-1)*N + i) = y(k) + b(k)*(s-x(k)) + c(k)*(s-x(k))^2 ...
                                 + d(k)*(s-x(k))^3;
        profile_Rd1((k-1)*N + i) = b(k) + 2*c(k)*(s-x(k)) ...
                                 + 3*d(k)*(s-x(k))^2;
        profile_Rd2((k-1)*N + i) = 2*c(k) ...
                                 + 6*d(k)*(s-x(k));
    end
end
profile_position(n*N+1) = x(end);
profile_R(n*N+1) = y(end);
profile_Rd1(n*N+1) = dy1;
profile_Rd2(n*N+1) = ddy1;

%% 三次样条差值结果展示
figure(1)
plot(profile_position,profile_R,'-',x,y,'*')
title('v^2-s 图')
figure(2)
plot(profile_position,profile_Rd1,'-',x,2*a,'*')
title('2a-s 图')
figure(3)
plot(profile_position,profile_Rd2,'-')
title('jerk/v-s图')

%% 拟合s-t图
profile_vel = sqrt(profile_R); 
t(1) = 0;
for k=2:N*n+1
    t(k) = t(k-1) + 2*ds/(profile_vel(k-1)+profile_vel(k));%路程除以速度等于时间
end

for k = 1:n+1
    t_0(k) = t((k-1)*N+1); %记录原始点对应的时间
end

dt = 0.01;
tt = t(1):dt:t(end);
sp = spline(t,profile_position,tt);
spd1 = diff(sp)/dt;
spd2 = diff(spd1)/dt;
spd3 = diff(spd2)/dt;
profile_acc = diff(profile_vel)./diff(t);
subplot(2,2,1)
plot(tt,sp,'-',t_0,x,'*');
title('s-t图')
subplot(2,2,2)
plot(tt(2:end),spd1,'-',t_0,v,'*')
title('v-t图')
subplot(2,2,3)
plot(tt(3:end),spd2,'-',t_0,a,'*')
title('a-t图')
subplot(2,2,4)
plot(tt(4:end),spd3)
title('jerk-t图')


