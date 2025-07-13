function [J_matrix] = jacobian_evalue(r_dot,r_r,r_sv)
%估算f相对于r_r的雅可比矩阵
%   此处显示详细说明
%% 估计法
% delta=1e3;
% deltaI=eye(3).*delta;
% r_r_delta=[r_r r_r r_r]+deltaI;
% f_value=f(r_dot,r_r,r_sv);
% f_j_1=(f(r_dot,r_r_delta(:,1),r_sv)-f_value)/delta;
% f_j_2=(f(r_dot,r_r_delta(:,2),r_sv)-f_value)/delta;
% f_j_3=(f(r_dot,r_r_delta(:,3),r_sv)-f_value)/delta;
% 
% J_matrix2=[f_j_1 f_j_2 f_j_3];
%% 精确法
sigma2=norm(r_r-r_sv);
sigma1=sigma2^3;
sigma3=r_dot.'*(r_r-r_sv);

jacobi=@(r_dot,r_r,r_sv) r_dot/sigma2-((r_r-r_sv)*sigma3/sigma1);
J_matrix=arrayfun(jacobi,r_dot,r_r,r_sv)';

end