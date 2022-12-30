clear all;

fxy = @(x,y) x+(3/2)*y.^2;
% marginal pdf
fx = @(x)integral(@(y)fxy,0,1);
fy = @(y)integral(@(x)fxy,0,1);

%calculate E[X]
fx1= @(x) x.^2.+(x./2);
Ex= integral(fx1,0,1);

%calculate E[X^2]
fx2= @(x) x.*(x.^2.+(x./2));
Ex2= integral(fx2,0,1);

%calculate var[X]
Var_X=Ex2-(Ex^2);

%calculate E[Y]
fy1= @(y) y./2 +(3/2)*y.^3;
Ey= integral(fy1,0,1);

%calculate E[Y^2]
fy2= @(y) y.*(y./2 +(3/2)*y.^3);
Ey2= integral(fy2,0,1);

%calculate var[Y]
Var_y=Ey2-(Ey^2);

%calculate E[XY]
fxy1=@(x,y) x.^2.*y+(3/2)*x.*y.^3;
Exy= integral2(fxy1,0,1,0,1);
cov_xy= Exy-(Ex*Ey);

%correlation matrix of U
R_U=[Ex2 Exy; Exy Ey2];
%covariance matrix of U
cov_U=[Var_X cov_xy; cov_xy Var_y];
disp(cov_U);
%cholesky factorization and Xs generation
Upper_tri=chol(cov_U,"upper");
Lower_tri=chol(cov_U,"lower");
Xs=Upper_tri'*randn(2,10000);
Cov_Xs=(cov(Xs'));
disp(Xs);

diff=cov_U-Cov_Xs;


