function [zmat] = zcalc(dopplershift,carrierfrequency,r_r,r_dot,r_sv,nSamples)
%  计算z,dopplershift:多普勒频移，carrierfrequency:载波
%   输入：r_dot：三维速度,3x1
%    r_init,位置初始值,3x1
%    r_sv，位置，3x1

zmat=zeros(nSamples,1);
for m=1:nSamples
    zmat(m)=dopplershift(m)*(3e8/carrierfrequency)-f(r_dot(:,m),r_r,r_sv(:,m));
end

end

