function [Area,Cx,Cy] = Centroid(x1,y1)
I = zeros(3,3);
Area = 0;
Cx = 0;
Cy = 0;
den = 0;

x = x1;
y = y1;

for i=1:1:length(x1)-1
    Area = Area + 1/2*(x(i)*y(i+1) - x(i+1)*y(i));
end
Area = abs(Area);

for i=1:1:length(x1)-1
    Cx = Cx + 1/(6*Area)*((x(i)+x(i+1))*(x(i)*y(i+1) - x(i+1)*y(i)));
    Cy = Cy + 1/(6*Area)*((y(i)+y(i+1))*(x(i)*y(i+1) - x(i+1)*y(i)));

end


end

