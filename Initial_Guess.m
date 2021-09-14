function [pos,vxp,vyp,vwp,alpha,delta,gamma,midx,midy] = Initial_Guess(pos,Vb,alpha,delta,gamma,zx,zy,xi,yi,time)
% Initial guess of values
global interval
global k

if k > 15
    [yf] = Poly_Expansion(time(k-3:k),Vb(k-3:k-1,:));
    pos(k,1) = pos(k-1,1) + yf(1)*interval;
    pos(k,2) = pos(k-1,2) + yf(2)*interval;
    pos(k,3) = pos(k-1,3) + yf(3)*interval;
    vxp = yf(1);
    vyp = yf(2);
    vwp = yf(3);
else
    pos(k,1) = pos(k-1,1) + Vb(k-1,1)*interval;
    pos(k,2) = pos(k-1,2) + Vb(k-1,2)*interval;
    pos(k,3) = pos(k-1,3) + Vb(k-1,3)*interval;
    vxp = Vb(k-1,1);
    vyp = Vb(k-1,2);
    vwp = Vb(k-1,3);
end


if k > 2 
    alpha(k) = alpha(k-1) + alpha(k-1)-alpha(k-2);
    delta(k) = delta(k-1) + delta(k-1)-delta(k-2);
    gamma(k) = gamma(k-1) + gamma(k-1)-gamma(k-2);
else
    delta(k) = 1/40;
    r = length(xi);
    pos(k,3) = 0.1*pi/180;
    for i=1:r
        zx(i,1) = [cos(pos(k,3)) -sin(pos(k,3))]*[xi(i);yi(i)] + pos(k,1);
        zy(i,1) = [sin(pos(k,3)) cos(pos(k,3))]*[xi(i);yi(i)] + pos(k,2);
    end
end

for i=1:length(zx)-1;
    midx(i,1) = (zx(i+1)+zx(i))/2;
    midy(i,1) = (zy(i+1)+zy(i))/2;
    del(i,1) = sqrt((zx(i)-zx(i+1))^2+(zy(i)-zy(i+1))^2);
    ss(i,1) = (zy(i+1)-zy(i))/del(i);
    cs(i,1) = (zx(i+1)-zx(i))/del(i);
end  



end

