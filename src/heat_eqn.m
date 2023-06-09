% size of domain (make it odd to have a point at 0.5,0.5)
N=5;
% number of design variables
d=2*N + 2*(N-2);
% design varibles
D=ones(d,1);
% location of cost function
l = floor(N/2)*N+floor(N/2)+1;
% solve analysis problem
[q0,dRdq]=solve_laplace(N,D);
% find sensitivities by finite difference
tic
D0=D;
eps=1e-2;
for k=1:d
 D=D0;
 D(k)=D0(k)+eps;
 [q1,dRdq1] = solve_laplace(N,D);
 dLdD(k) = (q1(l)-q0(l))/eps;
end
toc
display('Finite Difference Sensitivity')
% find senstivites using adjoint
tic
dLdq = zeros(N*N,1);
dLdq(l)=1;
% assemble dRdD knowing the discrete equation
% dRdQ*q-rhs(D) = 0;
dRdD  = zeros(N*N,d);
% lower plane
i=1;
k=1;
for j=1:N
 m=(i-1)*N+j;
 dRdD(m,k)=-1;
 k=k+1;
end
% upper plane
i=N;
for j=1:N
 m=(i-1)*N+j;
 dRdD(m,k)=-1;
 k=k+1;
end
% left plane
j=1;
for i=2:N-1
 m=(i-1)*N+j;
 dRdD(m,k)=-1;
 k=k+1;
end
% right plane
j=N;
for i=2:N-1
 m=(i-1)*N+j;
 dRdD(m,k)=-1;
 k=k+1;
end
% compute sensitivity using adjoint method
dLdD_adj = -dRdD'*inv(dRdq')*dLdq;
toc
display('Adjoint Sensitivity with full linearization')
% finite difference to find dRdD
tic
for k=1:d
 D=D0;
 D(k)=D0(k)+eps;
 rhs=form_rhs(N,D);
 dRdD_1(:,k)=(dRdq*q0-rhs)/eps;
end
dLdD_adj1 = -dRdD_1'*inv(dRdq')*dLdq;
toc
display('Adjoint Sensitivity with partial finite difference')
display('Norm of differences')
% check error of adjoint based sensitivities with finite difference
norm(dLdD_adj'-dLdD)
norm(dLdD_adj'-dLdD_adj1')
