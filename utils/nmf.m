function [W, H] = nmf(X, r)
% X - input matrix
% r - rank of embedding matrix
% max_iter 

max_iter=500;
% initialize
W = rand(size(X, 1), r);
H = rand(r, size(X, 2));

% iterate
for i = 1:max_iter
    % update H
    H = H .* (W.' * X) ./ (W.' * W * H);
    
    % update W
    W = W .* (X * H.') ./ (W * (H * H.'));
    
    W(W < 0) = 0;
    H(H < 0) = 0;
end
end