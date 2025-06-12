function M = matrix_vector_multiplication(A)
% A = [A1 A2 A3]
%     [A4 A5 A6]
%     [A7 A8 A9]
% M * B = A*B

nv = size(A,1);
n = sqrt(size(A,2));

I = repmat((1:nv)', [n,n]) + repelem((0:n-1)*nv, n*nv,1);
J = repmat((1:nv)', [n,n]) + repelem(repmat((0:n-1)'*nv, [1,n]), nv,1);
S = reshape(A, [n*nv,n]);
M = sparse(I, J, S, n*nv, n*nv);