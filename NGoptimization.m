function [r_r_est, a, err_hist] = NGoptimization(Z, f_d, r_sv, dot_r_sv, r_r_0, epsilon, max_iter)
% Newton-Gauss法求解多普勒定位最小二乘
% Z: 1x6 cell，每颗卫星观测
% f_d: 1x6 cell，每颗卫星f_d_i(k)
% r_sv: 1x6 cell，卫星位置
% dot_r_sv: 1x6 cell，卫星速度
% r_r_0: 3x1，初始坐标
% epsilon: 收敛阈值
% max_iter: 最大迭代次数

N_sv = numel(Z);
N_i = cellfun(@numel, Z);
N_tot = sum(N_i);

r_r = r_r_0;
a = zeros(N_sv,1);
err_hist = [];

for m = 1:max_iter
    z = zeros(N_tot,1);
    J = zeros(N_tot, 3+N_sv);
    cnt = 0;
    for i = 1:N_sv
        n = N_i(i);
        for k = 1:n
            cnt = cnt + 1;
            Zik = Z{i}(k);
            f_dik = f_d{i}(k);
            r_svik = r_sv{i}(:,k);
            dot_r_svik = dot_r_sv{i}(:,k);

            diff_r = r_r - r_svik;
            norm_diff_r = norm(diff_r);
            f_r = dot_r_svik' * diff_r / norm_diff_r;

            z(cnt) = Zik - f_r - f_dik*a(i);

            df_dr = (dot_r_svik' - f_r * (diff_r'/norm_diff_r)) / norm_diff_r;
            J(cnt, 1:3) = -df_dr;
            J(cnt, 3+i) = -f_dik;
        end
    end

    x_hat = (J'*J)\(J'*z);
    delta_r_r = x_hat(1:3);
    delta_a = x_hat(4:end);

    r_r = r_r + delta_r_r;
    a = a + delta_a;

    nu_z = z - J*x_hat;
    err = norm(nu_z)/norm(z);
    err_hist = [err_hist; err];

    if norm(delta_r_r) < epsilon || err < epsilon
        break
    end
end

r_r_est = r_r;
end

