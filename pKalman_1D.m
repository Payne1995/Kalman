clc;clear;close all
%% Kalman_1D
%   S 真实信号值 S[n] = A*S[n-1] + B*u[n-1] + w[n-1];
%   X 先验估计  X[n] = A*Y[n-1] + B*u[n-1];
%   Z 测试值   Z[n] = H*S[n] + v[n];
%   Y 后验估计  Y[n] = X[n] - K[n](Z[n] - X[n]);
%   P1 先验误差协方差   P1 = A*P0*A' +  R;
%   K 卡尔曼增益     K = (P1*H')/(H*P1*H' + Q)
%   P0 后验误差协方差   P0 = (I - K*H)P1;
%   过程噪声 R
%   测量噪声 Q
%% 一维卡尔曼赋值
%   A = 0.5；
%   B = 2;
%   H = 3;

%% Init
N = 1000;
A = 0.5;
B = 0.5;
H = 1;
S = zeros(1,N);
Z = zeros(1,N);
u = linspace(1,100,N);
R = 0.1;            % 过程噪声
Q = 10;             % 测量噪声
w = sqrt(R).*randn(N,1);
v = sqrt(Q).*randn(N,1);

X = zeros(1,N);
Y = zeros(1,N);
%% 信号与测量值
for i = 2:N
   S(i) = A * S(i-1) + B * u(i-1) + w(i-1); 
   Z(i) = H * S(i) + v(i); 
end
%% 最优估计
Y(1) = S(1);
P0 = eye(1);
for i = 2:N
    Xn = A * Y(i-1) + B * u(i-1);
    P1 = A * P0 * A' + R;
    K = (P1 * H') * inv(H * P1 * H' + Q);
    Y(i) = Xn + K*(Z(i) - H*Xn);
    P0 = (eye(1) - K*H)*P1;
end

%% plot
figure(1)
hold on
plot(Z(2:N),'b')
plot(S(2:N),'r')
legend('Measured ','Actual');
xlabel('time, n')
ylabel('Signal')
grid on
figure(2)
hold on
plot(Y(2:N),'b')
plot(S(2:N),'r')
legend('Estimated ','Actual');
xlabel('time, n')
ylabel('Signal')
grid on