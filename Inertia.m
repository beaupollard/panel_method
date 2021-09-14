function [I] = Inertia(x1,y1,Cx,Cy)

global mass

I = zeros(3,3);
J = 0;
x = x1;
y = y1;
den = 0;
x = x - Cx;
y = y - Cy;
for i=1:1:length(x1)-1
    xp = x(i);
    xp1 = x(i+1);
    yp = y(i);
    yp1 = y(i+1);
    
    a = xp*yp1 - xp1*yp;
    I(1,1) = mass;
    I(2,2) = mass;
    I(3,3) = I(3,3) + mass/6*(xp*yp1-yp*xp1)*(xp^2+yp^2+xp*xp1+yp*yp1+xp1^2+yp1^2);
    den = den+(xp*yp1-yp*xp1);
end
I(3,3) = I(3,3)/den;

end

