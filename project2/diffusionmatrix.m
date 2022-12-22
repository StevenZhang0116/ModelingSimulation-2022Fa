function [A, B, G] = diffusionmatrix(dt, R, N, D, L, I, C)
% Calculating diffusion matrix for axonal propagation 

% Input 
%  dt  : size of time step [ms] (scalar)
%  R : Axial resistivity [MOhm.um] = [Ohm.m] (scalar)
%  N : Number of nodes (scalar)
%  D : Axonal diameter [um] (vector)
%  L : Nodal length [um] (vector)
%  I : Internodal length [um] (vector)(use I=0 for unmyelinated axons)
%  C : Nodal capacitance [uF] (vector)
%
% Output
%  A, B : diffusion matrices [ A*u(t+1) = B*u(t) ] 
%  G : axial resistance [mS]

% calculating relevant values 
% make full output matrix
Af = zeros(N,N);
Bf = zeros(N,N);

% Axonal resistance
R0 = R * (L(1:end-1)/2) ./ ( pi * (D(1:end-1)/2).^2 ); %[MOhm]
R1 = R *  I             ./ ( pi * (D(1:end-1)/2).^2 ); %[MOhm]
R2 = R * (L(2:end)  /2) ./ ( pi * (D(2:end)  /2).^2 ); %[MOhm]
g = 1./(R0+R1+R2)*1e-3; % [mS] axial conductance; length(g)=N-1

% assigning data values to diffusion matrix
% one end 
r = g(1) * dt / C(1) / 2; % [mS*ms/uF]=[no unit]
Af(1,1) = 1+r;
Af(1,2) =  -r;
Bf(1,1) = 1-r;
Bf(1,2) =  +r; 
% intermediate nodes 
for i = 2:N-1
 r1 = g(i-1) * dt / C(i) / 2; 
 r2 = g(i)   * dt / C(i) / 2; 
 Af(i,i-1) = -r1;
 Af(i,i)   = 1+r1+r2;
 Af(i,i+1) = -r2;
 Bf(i,i-1) = +r1;
 Bf(i,i)   = 1-r1-r2;
 Bf(i,i+1) = +r2;
end
% other end 
r = g(N-1) * dt / C(N) / 2;
Af(N,N-1) =  -r;
Af(N,N)   = 1+r;
Bf(N,N-1) =  +r; 
Bf(N,N)   = 1-r;

% converting into sparse matrix (to reduce time for computing diffusion) 
A = sparse(Af);
B = sparse(Bf);
G = g; 

