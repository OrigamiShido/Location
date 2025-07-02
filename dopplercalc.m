function [frequencyRate,position,velosity]=dopplercalc(targetTime)
%% 定义数值
% tic
disp('begin...')
latitude=30.5288888;
longitude=114.3530555;
altitude=56;
minelevation=0;
durationtimeSeconds=120;
starttime=targetTime;
sampletime=1;

%% 模拟环境
disp('creating satellitescenario object...')
% 创建图窗
sc = satelliteScenario(starttime,starttime+seconds(durationtimeSeconds),sampletime);

%% 地面站
disp('creating groundstation...')
% 创建地面站
gs=groundStation(sc,Name='WHU',Latitude=latitude,Longitude=longitude,Altitude=altitude,MinElevationAngle=minelevation);

%% 导入卫星
disp('importing satellites...')
% 创建和读取卫星，渲染轨道
sat=satellite(sc,[pwd,'\gp2.tle'],OrbitPropagator="sgp4");
%% 预报多普勒频移
disp('Predicting Doppler shift...')
carrierFrequency=11.325e9;
[frequencyShift,timeOut,dopplerInfo] = dopplershift(sat,gs,Frequency=carrierFrequency);
frequencyRate = dopplerInfo.DopplerRate;
% relativeVelocity = dopplerInfo.RelativeVelocity;

%% 做表，删除nan行
rowname=string(starttime:seconds(sampletime):starttime+seconds(durationtimeSeconds-1));
frequencyRate=array2table(frequencyRate','RowNames',rowname,'VariableNames',sat.Name);
frequencyRate(:,all(isnan(frequencyRate{:,:})))=[];

%% 获取位置
disp('getting positions and predicting...')

[position,velosity]=states(sat,CoordinateFrame='ecef');
% [position_lla,velocity_lla]=states(sat,CoordinateFrame='geographic');
%% 预报

% 获取地面站 ECEF 坐标
% gsLLA = [latitude, longitude, altitude];
% gsECEF = lla2ecef(gsLLA);
% % 转换为表
% rowname=string(starttime:seconds(sampletime):starttime+seconds(durationtimeSeconds));

end
