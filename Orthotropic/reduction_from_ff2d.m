function [k21,Reduction] = reduction_from_ff2d(Src, param, ang, omega, Edge_jump, v2t)

% Compute jumps per edge from cross field
k21 = ones(Src.ne,1);
[~,k21i] = min(abs(exp(1i*ang(param.E2T(param.ide_int,2))+(0:3)*1i*pi/2+1i*(omega(param.ide_int)-param.para_trans(param.ide_int))) - exp(1i*ang(param.E2T(param.ide_int,1)))), [], 2);
k21(param.ide_int) = k21i;

% Accumulate jump around vertices
k21T = mod(reshape(Edge_jump*(k21-1), Src.nf,[]), 4);

% Sign bits
s = (-1).^k21T;

% Build reduction matrix
v2t_smooth = sparse(reshape((1:3*Src.nf)', [Src.nf,3]), Src.T, 1, 3*Src.nf, Src.nv);
Reduction = blkdiag(v2t_smooth, spdiags(s(:), 0, 3*Src.nf, 3*Src.nf)*v2t);
