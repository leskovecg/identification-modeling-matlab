function z=xcorrp(x,y)
% function z=xcorrp(x,y)
%
% Function calculates cross correlation of periodic functions.
% Calculation is done only for shifts 0, 1, 2, ... length(u)-1.
% The result is periodic, and the correlation for neagtive shifts can be obtained from posiive ones.
%
% x and y must be (column or row) vectors of the same dimensions.
% z is of the same dimension as x and y.

[n1,n2] = size(x);
[n3,n4] = size(y);

if (n1 ~= 1) && (n2 ~= 1)
    error('The first argument must be (column or row) vector.')
end

if (n3 ~= 1) && (n4 ~= 1)
    error('The second argument must be (column or row) vector.')
end

if (n1 ~= n3) || (n2 ~= n4)
    error('The dimensions of the input vectors should be the same.')
end

N = length(x);
z = NaN(size(x)); % initialization of the result with not-a-number

x2 = x(:); % auxiliary column vector for x - operator (:) makes column vector
y2 = [y(:); y(:)]; % auxiliary column vector for two periods of  y

for i=1:N
    z(i)=x2' * y2(i:i+N-1); % x is fixed, y is shifted
end

z = z/N; % normalization with the length
