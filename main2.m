%% 
% 示例表格
T = table([1;2;3], [4;5;6], 'VariableNames', {'A', 'B'}, 'RowNames', {'X', 'Y', 'Z'});

% 已知行名和列名
rowName = 'Y';
colName = 'B';

% 找索引
rowIdx = find(strcmp(T.Properties.RowNames, rowName));
colIdx = find(strcmp(T.Properties.VariableNames, colName));

% 数值索引
disp(['行索引: ', num2str(rowIdx), ', 列索引: ', num2str(colIdx)]);
% 取出该点的值
value = T{rowIdx, colIdx};
disp(['对应数值: ', num2str(value)]);

%% 

% Step 2: 整理数据
sv_names = frequencyRate.Properties.VariableNames; % 卫星名称
time_names = frequencyRate.Properties.RowNames;    % 时间字符串
N_sv = numel(sv_names);
N_obs = numel(time_names);


% 假设position和velosity是 3 x N_sv x N_obs 或 3 x N_obs x N_sv 格式
% 这里假定为 3 x N_obs x N_sv
Z = cell(1, N_sv);
f_d = cell(1, N_sv);
r_sv = cell(1, N_sv);
dot_r_sv = cell(1, N_sv);

for i = 1:N_sv
    % 多普勒频移观测，N_obs x 1
    Z{i} = frequencyRate{:, sv_names{i}};
    % f_d 全1
    f_d{i} = ones(N_obs, 1);
    % 卫星位置和速度，3 x N_obs
    r_sv{i} = squeeze(position(:, :, i));    % 3 x N_obs
    dot_r_sv{i} = squeeze(velosity(:, :, i));% 3 x N_obs
end

% Step 3: 随机生成初始ECEF坐标
% 以地心为中心，半径在[6371km, 6371km+1000km]范围内随机一个点为初值（可按实际需求调整）
earth_radius = 6371000;
r_r_0 = earth_radius * (rand(3,1)-0.5)*2 + [0;0;0];

% Step 4: 设定收敛参数
epsilon = 1e-6;
max_iter = 20;

% Step 5: 调用最小二乘法主函数
[r_r_est, a_est, err_hist] = NGoptimization(Z, f_d, r_sv, dot_r_sv, r_r_0, epsilon, max_iter);

% Step 6: 显示结果
disp("估计的接收机坐标(ECEF):");
disp(r_r_est);
disp("估计的偏差参数:");
disp(a_est);
disp("误差收敛历史:");
disp(err_hist);


