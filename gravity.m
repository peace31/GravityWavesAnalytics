function [u2,h2,eta2] = gravity(u1,h1,eta1,N,dx,dt,g,H);
% USAGE: [u2,h2,eta2] = gravity(u1,h1,eta1,N,dx,dt,g,H);
%
% J. S. Sabato - April, 2018
%
% Inputs:
% u1            - flow speed at nodes, time level 1 (N+1) X 1
% eta1          - surface elevation at nodes, time level 1 (N+1) X 1
% h1            - surface elevation at nodes, time level 1 (N+1) X 1
% N,dx,dt,g,H   - array size, grid size, time step, gravity, 
%                 mean depth 1 X 1
% 
% Outputs:
% u2            - flow speed at nodes, time level 2 after dt  (N+1) X 1
% eta2          - surface elevation at nodes, time level 2 after dt (N+1) X 1
% h2            - surface elevation at nodes, time level 2 after dt (N+1) X 1
%

% check sizes of input arguments, return error

if sum([size(dx)~=[1 1] size(N)~=[1 1] size(dt)~=[1 1] size(g)~=[1 1] size(H)~=[1 1]])
    error('ERROR: one or more of dx, N, dt, g, H is not scalar');
elseif sum([size(u1)~=[N+1 1] size(eta1)~=[N+1 1] size(h1)~=[N+1 1]])
    error('One or more vector input arguments are wrong size');    
else 
end 

% construct A,b 
A = zeros(2*(N+1));
b = zeros(2*(N+1),1);
    % u on boundaries
    b(1) = 0;
    A(1,1) = 1;
    b(N+1) = 0;
    A(N+1,N+1) = 1;
    % eta on boundaries
    b((N+1)+1) = 0;
    A((N+1)+1,(N+1)+1) = 1;
    A((N+1)+1,(N+1)+2) = -1;
    b(2*(N+1)) = 0;
    A(2*(N+1),2*(N+1)) = 1;
    A(2*(N+1),2*(N+1)-1) = -1;
for n=2:N
    % u interior
    b(n) = u1(n);
    A(n,n) = 1;
    A(n,(N+1)+(n-1)) = -g*dt/(2*dx);
    A(n,(N+1)+(n+1)) = g*dt/(2*dx);    
    % eta interior
    b((N+1)+n) = eta1(n);
    A((N+1)+n,(N+1)+n) = 1;
    A((N+1)+n,n-1) = -H*dt/(2*dx);
    A((N+1)+n,n+1) = H*dt/(2*dx);  
end
% solve Ax=b
x = A\b;
% decompose u, eta
u2 = x(1:N+1);
eta2 = x(N+2:2*(N+1));
h2 = eta2+H;

end