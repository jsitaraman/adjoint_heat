% perform laplace equation on a uniform Cartesian
% grid with prescribed boundary conditions given in
% the vector D
%
% Jay Sitaraman
% 6/7/2023
function rhs=form_rhs(N,D)
% Solve [dRdq] [q] = rhs(D)
% direct inversion
rhs=zeros(N*N,1);
%
% set RHS to boundary conditions
% 
k=1;
j=1;
for i=1:N
  m=(i-1)*N+j;
  rhs(m)=D(k);
  k=k+1;
end
%
j=N;
for i=1:N
 m=(i-1)*N+j;
 rhs(m)=D(k);
 k=k+1;
end
%
i=1;
for j=2:N-1
 m=(i-1)*N+j;
 rhs(m)=D(k);
 k=k+1;
end
%
i=N;
for j=2:N-1
 m=(i-1)*N+j;
 rhs(m)=D(k);
 k=k+1;
end
