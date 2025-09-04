function C = OuterProduct(A, B)  % version 5
C = reshape(A(:) * B(:).', [size(A), size(B)]);
end