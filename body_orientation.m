function [n1,t1,del,midx,midy,zx,zy] = body_orientation(pos,xi,yi)

global k

r = length(xi);

if k == 2
    pos(k,3) = 1*0.01;
end
for i=1:r
    zx(i,1) = [cos(pos(k,3)) -sin(pos(k,3))]*[xi(i);yi(i)] + pos(k,1);
    zy(i,1) = [sin(pos(k,3)) cos(pos(k,3))]*[xi(i);yi(i)] + pos(k,2);
end

for i=1:length(zx)-1
    midx(i,1) = (zx(i+1)+zx(i))/2;
    midy(i,1) = (zy(i+1)+zy(i))/2;
    del(i,1) = sqrt((zx(i)-zx(i+1))^2+(zy(i)-zy(i+1))^2);
    ss(i,1) = (zy(i+1)-zy(i))/del(i);
    cs(i,1) = (zx(i+1)-zx(i))/del(i);
end

n1 = [-ss(:) cs(:)];
t1 = [cs(:) ss(:)];   


end

