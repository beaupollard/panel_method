function [yf] = Poly_Expansion(x,y)

[r,c] = size(y);
xf = x(end);

for h=1:c
    sum = 0;
    for i = 1:r
        prod = 1;
        for j = 1:r
            if j ~= i
                prod = prod*(xf-x(j))/(x(i)-x(j));
            end
        end
        sum = sum+prod*y(i,h);
    end
    yf(h) = sum;
end

