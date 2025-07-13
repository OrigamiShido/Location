function r_r=newton_gauss(r_dot,r_r,r_sv,dopplerShift,carrierFrequency,epsilon,max_iter)
% 牛顿-高斯法解定位方程：
% 输入：r_dot:卫星的ECEF速度，3x采样点数量
% r_r:初始位置，3x1
% r_sv:卫星的ECEF位置，3x采样点数目
% dopplerShift：多普勒频移，1x采样点数目
% carrierFrequency:常数，载波频率

numSample=size(r_dot,2);
nu_history=[];
% r_r_init=r_r;
% J=Jcalc(r_dot,r_r,r_sv,numSample);

for m=1:max_iter
    % 更新z,J

    J=Jcalc(r_dot,r_r,r_sv,numSample);
    z=zcalc(dopplerShift,carrierFrequency,r_r,r_dot,r_sv,numSample);

    % 计算x_hat

    x_hat=pinv((J')*J)*(J')*z;

    % r_r前进步长

    deltax=[x_hat(1) x_hat(2) x_hat(3)]';

    r_r=r_r+1000.*deltax;
    
    ecef2lla(r_r')

    % 计算nuz

    nu=z-J*x_hat;
    if m==1
        nu_1=nu;
    end
    err=norm(nu)/norm(nu_1);
    nu_history=[nu_history err];

    if (norm(deltax)<=epsilon||err<=epsilon)
        break;
    end
end
plot(nu_history);
end