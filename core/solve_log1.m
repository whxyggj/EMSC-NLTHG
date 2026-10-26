function [E] = solve_log1(W,lambda)
% log(1+|x|)
n = size(W,2);
E = W;
for i=1:n
    E(:,i) = solve_log11(W(:,i),lambda);
end
end

function [x] = solve_log11(w,lambda)

nw = norm(w);
a = (nw-1+((nw+1)^2-4*lambda)^0.5)/2;
if nw>lambda
    x = (1-lambda/(lambda+a+a^2))*w;
else
    x = zeros(length(w),1);
end
end
