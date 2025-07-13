%% 数值定义
targetTime=datetime(2025,7,7,16,45,0,'TimeZone',hours(8));
latitude_true=30.5288888;
longitude_true=114.3530555;
altitude_true=56;
carrierfrequency=11.325e9;
lambda=3e8/carrierfrequency;
%% 正过程：获得目标时间内的多普勒频移
[frequencyRate,position,velosity]=dopplercalc(targetTime);
%%
frequencyRef=frequencyRate;
frequencyRef(:,all(isnan(frequencyRef{:,:})))=[];

%%
% 已知行名和列名
valName = frequencyRef.Properties.VariableNames;

% 找索引
[~,valIdx] = ismember(valName,frequencyRate.Properties.VariableNames);

val=valIdx([1 5 6 7 8 9]);

%% 反过程：构建方程并求解

% 选取值
dopplerShift=frequencyRate{1,val}';
r_sv=squeeze(position(:,1,val));
r_dot=squeeze(velosity(:,1,val));

% 假定真值坐标
r_r_true=lla2ecef([latitude_true longitude_true altitude_true]);


% 地球半径（近似，单位：米）
R = 6378137;

% 2. 随机方位角和距离
theta = rand() * 2 * pi;             % 方位角 [0, 2π]
phi = acos(1 - rand() * (1 - cos(100000/R)));  % 球面距离角度，确保均匀分布在球帽范围

% 新点相对中心点的球面坐标（单位弧度）
d = 100000; % 最大距离 (米)
% 推荐更标准方式如下
delta_sigma = acos(1 - rand() * (1 - cos(d/R)));

% 3. 计算新点经纬度
lat0 = deg2rad(latitude_true);
lon0 = deg2rad(longitude_true);

% 新点纬度
lat_new = asin(sin(lat0)*cos(delta_sigma) + cos(lat0)*sin(delta_sigma)*cos(theta));

% 新点经度
lon_new = lon0 + atan2(sin(theta)*sin(delta_sigma)*cos(lat0), cos(delta_sigma) - sin(lat0)*sin(lat_new));

% 新点高度（与原点一样）
alt_new = altitude_true;

% 4. 转回ECEF
lla_new = [rad2deg(lat_new), rad2deg(lon_new), alt_new];
r_r_0 = lla2ecef(lla_new)';
 
r_r=newton_gauss(r_dot,r_r_0,r_sv,dopplerShift,carrierfrequency,1e-4,200000);

r_r_true_lla=ecef2lla(r_r_true);
r_r_0_lla=ecef2lla(r_r_0');

r_r_lla=ecef2lla(r_r');

[xn,ye,zup]=ecef2enu(r_r(1),r_r(2),r_r(3),latitude_true,longitude_true,altitude_true,wgs84Ellipsoid);

err_2d=norm([xn ye]);

err_3d=norm([xn,ye,zup]);

% plot(error_history)