function [value] = f(r_dot,r_r,r_sv)
%CALCF 计算f
%   输入：r_dot：三维速度,3x1
%    r_r,位置估计,3x1
%    r_sv，位置，3x1
value=(r_dot'*(r_r-r_sv))/norm(r_r-r_sv);

end

