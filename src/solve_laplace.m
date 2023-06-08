% perform laplace equation on a uniform Cartesian
% grid with prescribed boundary conditions given in
% the vector D
%
% Jay Sitaraman
% 6/7/2023
function [q,dRdq]=solve_laplace(N,D)
% Solve [dRdq] [q] = rhs(D)
% direct inversion
dRdq=zeros(N*N,N*N);
rhs=zeros(N*N,1);
%
% set RHS to boundary conditions
% and LHS to 1 corresponding to the
% row of the BC
%
k=1;
j=1;
for i=1:N
  m=(i-1)*N+j;
  rhs(m)=D(k);
  dRdq(m,m)=1;
  k=k+1;
end
%
j=N;
for i=1:N
 m=(i-1)*N+j;
 rhs(m)=D(k);
 dRdq(m,m)=1;
 k=k+1;
end
%
i=1;
for j=2:N-1
 m=(i-1)*N+j;
 rhs(m)=D(k);
 dRdq(m,m)=1;
 k=k+1;
end
%
i=N;
for j=2:N-1
 m=(i-1)*N+j;
 rhs(m)=D(k);
 dRdq(m,m)=1;
 k=k+1;
end
% setup the Jacobian matrix
for i=2:N-1
 for j=2:N-1
  ij=(i-1)*N+j;
  im1j=ij-N;
  ip1j=ij+N;
  ijm1=ij-1;
  ijp1=ij+1;
  dRdq(ij,ij)=1;
  dRdq(ij,im1j)=-0.25;
  dRdq(ij,ip1j)=-0.25;
  dRdq(ij,ijm1)=-0.25;
  dRdq(ij,ijp1)=-0.25;
 end
end
% solve laplace by direct inversion
q=inv(dRdq)*rhs;
