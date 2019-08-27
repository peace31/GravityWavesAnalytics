function [u2,h2,eta2,phi2] = nonlinear(u1,h1,eta1,phi1,N,dx,dt,g,H);
% USAGE: [u2,h2,eta2,phi2] = nonlinear(u1,h1,eta1,phi1,N,dx,dt,g,H)
%
% J.S. Sabato, April 2018
%
% Inputs:
% u1            - old flow speed (N+1) X 1
% h1            - old fluid depth (N+1) X 1
% eta1          - old fluid depth (N+1) X 1
% phi1          - old tracer concentration (N+1) X 1
% N             - array size 1 X 1
% dx            - cell width / distance between nodes 1 X 1
% dt            - time step 1 X 1
% g, H          - gravity acceleration and mean fluid depth 1 X 1
% 
% Outputs:
% u2            - new flow speed (N+1) X 1
% h2            - new fluid depth (N+1) X 1
% eta2          - new fluid depth (N+1) X 1
% phi2          - new tracer concentration (N+1) X 1

	[u1,h1,eta1] = gravity(u1,h1,eta1,N,dx,dt,g,H);

    % cell-center values from cell-edge values
    % flow and surface
    uc1 = (u1(1:N)+u1(2:N+1))/2;
    hc1 = (h1(1:N)+h1(2:N+1))/2;
    % mass flux
    U1 = u1.*h1;
    Uc1 = (U1(1:N)+U1(2:N+1))/2;    
    % tracer
    phic1 = (phi1(1:N)+phi1(2:N+1))/2;
        
    %
    %
    % advect mass, mass flux and tracer
    [hc2] = FCT(hc1,uc1,dx,N,dt);
    [Uc2] = FCT(Uc1,uc1,dx,N,dt);
    [phic2] = FCT(phic1,uc1,dx,N,dt);
    uc2 = Uc2./hc2;
    %
    %
    %
        
    % edge-values from cell-centered values
    % surface
    h2(2:N) = (hc2(1:N-1)+hc2(2:N))/2;
    % BCs
    h2(1) = h2(2);
    h2(N+1) = h2(N);
    % momentum
    u2(2:N) = (uc2(1:N-1)+uc2(2:N))/2;
    % BCs
    u2(1) = 0;
    u2(N+1) = 0;
    % tracer
    phi2(2:N) = (phic2(1:N-1)+phic2(2:N))/2; 
    % BCs
    phi2(1) = phi2(2);
    phi2(N+1) = phi2(N); 
    
    eta2 = h2-H;
end

function [q2]=FCT(q1,u,dx,N,dt);
% USAGE: [q2]=FCT(q1,u,dx,N,dt)
%
% J.S. Sabato, April 2018
%
% 1-D advection of tracer following
% the flux conservation law dq/dt + d(uq)/dx = 0
% by flux corrected transport/flux limiter method
% with superbee limiter, no-flux boundary conditions
%
% q1: initial distribution of q (cell-centered)
% q2: distribution after dt seconds (cell-centered)
% u:  advecting flow speed (cell-centered)
% dt: time step
% dx: cell width
% N:  length of arrays
% N cells, N+1 edges/fluxes
% umh = flow speed at cell edge n-1/2
%
    f(1) = 0; 
for n=2:N
    if n==N
        ii = [N-2 N-1 N N];
        umh = (u(N-1)+u(N))/2;
    elseif n==2
        ii = [1 1 2 3];
        umh = (u(1)+u(2))/2;
    else
        ii = [n-2 n-1 n n+1];
        umh = (u(n-1)+u(n))/2;
    end
    [f(n)] = flux(umh,q1(ii),dx,dt);
end
    f(N+1) = 0;   
q2 = q1 - (dt/dx).*(f(2:N+1)'-f(1:N)');
end

function [f] = flux(u,q,dx,dt);
% USAGE: [f] = flux(u,q,dx,dt)
%
% J.S. Sabato, April 2018
%
% u = left-cell-edge flow u_{n-1/2}
% q = flux-conserved quantity (4-point stencil, cell centered)
% dx, dt = cell width, time interval
% superbee flux limiter
    if u>=0
        r = (q(2)-q(1))/(q(3)-q(2));
    else
        r = (q(4)-q(3))/(q(3)-q(2));
    end   
    phi  = max([0 min(1,2*r) min(2,r)]);
    f  =     (                                              ...
                     0.5*(  (1+sign(u))*q(2)*u              ... 
                           +(1-sign(u))*q(3)*u  )           ... 
             + 0.5*abs(u)*(1-abs(u*dt/dx))*phi*(q(3)-q(2))  ...
             ); 
end