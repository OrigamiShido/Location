z=zeros(6,1);
x=zeros(9,1);
J=zeros(6,7);
carrierFrequency=1.525e9;

epsilon=1e-6;
max_iter=20;
for m=1:max_iter
% 更新z,J

J=J(r_dot,r_r,r_sv);
z=z(dopplerShift,carrierFrequency,r_r,r_dot,r_sv);

% 计算x_hat

x_hat=(J.'*J)^-1*J.'*z;

% r_r前进步长

deltax=[x_hat(1) x_hat(2) x_hat(3)];

r_r=r_r+deltax;

% 计算nuz

nu=z-J*x;

    if (norm(deltax)<epsilon)
        break;
    end
end

% 计算误差
