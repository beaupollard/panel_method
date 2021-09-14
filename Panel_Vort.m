function [Ax,Ay,Bx,By] = Panel_Vort(xi,yi,xf,yf)
%Affect of a panel on a vortex
%% xi, yi are the points being influenced and xf, yf are the points influencing
ri = length(xi);
rf = length(xf);

for i=1:1:rf-1
    midx(i) = (xf(i+1)+xf(i))/2;
    midy(i) = (yf(i+1)+yf(i))/2;
    Dsf(i) = sqrt((xf(i+1)-xf(i))^2+(yf(i+1)-yf(i))^2);
    ssf(i) = (yf(i+1)-yf(i))/Dsf(i);
    csf(i) = (xf(i+1)-xf(i))/Dsf(i);
end

Ax = zeros(ri,rf-1);
Ay = zeros(ri,rf-1);   
Bx = zeros(ri,rf-1);
By = zeros(ri,rf-1);   

for i=1:1:rf-1
    xp = [];
    yp = [];

    for j=1:1:ri
        xp2 = xi(j) - midx(i);
        yp2 = yi(j) - midy(i);
        xp(j,1) = [csf(i) ssf(i)]*[xp2;yp2];
        yp(j,1) = [-ssf(i) csf(i)]*[xp2;yp2];
    end

    for j=1:1:ri
        Vt = 0;
        Vn = 0;
        Vtg = 0;
        Vng = 0;
        
        midxp = xp(j);
        midyp = yp(j);
        nor = 4*pi;
        num = (midxp + Dsf(i)/2)^2+midyp^2;
        den = (midxp - Dsf(i)/2)^2+midyp^2;
        Vt = log(num/den)/(nor);
        Vng = -Vt;
        
        num = midyp*Dsf(i);
        den = (midxp^2 + midyp^2 - (Dsf(i)/2)^2);
        Vn = 2*atan2(num,den)/(nor);
        Vtg = Vn;
        
        Ax(j,i) = Vt*csf(i) - Vn*ssf(i);
        Ay(j,i) = Vt*ssf(i) + Vn*csf(i);
        Bx(j,i) = Vtg*csf(i) - Vng*ssf(i);
        By(j,i) = Vtg*ssf(i) + Vng*csf(i);

    end
end

end

