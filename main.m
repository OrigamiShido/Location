%% 数值定义
targetTime=[datetime(2025,5,17,15,39,57,698,'TimeZone',hours(8)) datetime(2025,5,17,15,44,29,94,'TimeZone',hours(8))];
latitude_true=30.5288888;
longitude_true=114.3530555;
altitude_true=56;
carrierfrequency=11.325e9;
%% 正过程：获得目标时间内的多普勒频移
[frequencyRate,position,velosity]=dopplercalc(datetime(2025,5,17,15,39,57,698,'TimeZone',hours(8)));

%% 反过程：构建方程并求解

% 选取值
dopplershift=frequencyRate{1,[1 4 6 7 8 10]};
r_sv=squeeze(position(:,1,[1 4 6 7 8 10]));
r_dot=squeeze(velosity(:,1,[1 4 6 7 8 10]));

% 假定真值坐标
r_r_true=lla2ecef([latitude_true longitude_true altitude_true]);

% Step 3: 设置初始ECEF坐标（更合理的初始值）
% 在真实位置附近随机偏移，而不是完全随机
r_r_true_initial = lla2ecef([latitude_true longitude_true altitude_true]);
% 在真实位置周围1km范围内随机偏移
offset_range = 1000; % 1km偏移范围
r_r_0 = r_r_true_initial + (rand(3,1)-0.5)*2*offset_range;

r_r=newton_gauss(r_dot,r_r_0,r_sv,dopplershift,carrierfrequency,1e-6,2000);

r_r_lla=ecef2lla(r_r');

error=r_r_lla-[latitude_true longitude_true altitude_true];
error_norm=norm(error);

% 输出结果
fprintf('\n=== 定位结果 ===\n');
fprintf('真实位置: 纬度=%.6f°, 经度=%.6f°, 高度=%.1fm\n', latitude_true, longitude_true, altitude_true);
fprintf('估计位置: 纬度=%.6f°, 经度=%.6f°, 高度=%.1fm\n', r_r_lla(1), r_r_lla(2), r_r_lla(3));
fprintf('误差: 纬度=%.6f°, 经度=%.6f°, 高度=%.1fm\n', error(1), error(2), error(3));
fprintf('总误差范数: %.6f\n', error_norm);

% 计算水平距离误差（更直观）
R_earth = 6371000; % 地球半径
lat_error_m = error(1) * pi/180 * R_earth;
lon_error_m = error(2) * pi/180 * R_earth * cos(latitude_true * pi/180);
horizontal_error = sqrt(lat_error_m^2 + lon_error_m^2);
fprintf('水平位置误差: %.2f 米\n', horizontal_error);
fprintf('高度误差: %.2f 米\n', error(3));