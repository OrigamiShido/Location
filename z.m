function [value] = z(dopplershift,carrierfrequency,r_r,r_dot,r_sv)
%  计算z,dopplershift:多普勒频移，carrierfrequency:载波
%   输入：r_dot：三维速度,3x1
%    r_init,位置初始值,3x1
%    r_sv，位置，3x1

Z=dopplershift*(3e8/carrierfrequency);
value=Z-f(r_dot,r_r,r_sv);

end

