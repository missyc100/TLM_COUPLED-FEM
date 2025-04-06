function [U1,fl1] = solve(Keff,Feq)
A100=sparse(size(Keff,1),1,0);
%%   均衡—使得更加对角占优
[P_keff,R_keff,C_keff]=equilibrate(Keff);
B_keff=R_keff*P_keff*Keff*C_keff;
clear Keff
%%   排序—最小化矩阵因子中的非零元素

% tic
[ddd1]=symrcm(B_keff);% symrcm   dissect
B_Keff1=(B_keff(ddd1,ddd1));
% toc
clear B_keff

%% ilu分解，建立gmres预条件子
setup = struct('type','ilutp','droptol',1e-3,'thresh',0);
% tic
[L2,U2] = ilu(B_Keff1,setup);
% t1=toc
maxit=1e3;
tol=1e-8;
restart=4;
dd=R_keff*P_keff*Feq;
%  toc
%     [A1,fl1,rr1,it1,rv1] = gmres(Keff1,Feq(:,1),restart,tol,maxit,[L],[U],A0);%     A1=dKeff\Feq;
%     A100=luB_Keff1\dd(ddd1,:);
[A100,fl1,rr1,it1,rv1] = gmres(B_Keff1,dd(ddd1,:),restart,tol,maxit,[L2],[U2],A100);%     A1=dKeff\Feq;
A10(ddd1,1)=A100;
U1=C_keff*A10;
end


