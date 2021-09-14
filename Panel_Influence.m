function [An,At,Ax,Ay,Bn,Bt,Bx,By] = Panel_Influence(xi,yi,xf,yf)
%Influence from panels
%% xi, yi are the points being influenced and xf, yf are the points influencing
ri = length(xi);
rf = length(xf);

for i=1:1:ri-1
    Dsi(i) = sqrt((xi(i+1)-xi(i))^2+(yi(i+1)-yi(i))^2);
    ss(i) = (yi(i+1)-yi(i))/Dsi(i);
    cs(i) = (xi(i+1)-xi(i))/Dsi(i);
    midxi(i) = (xi(i) + xi(i+1))/2;
    midyi(i) = (yi(i) + yi(i+1))/2;
end

for i=1:1:rf-1
    midx(i) = (xf(i+1)+xf(i))/2;
    midy(i) = (yf(i+1)+yf(i))/2;
    Dsf(i) = sqrt((xf(i+1)-xf(i))^2+(yf(i+1)-yf(i))^2);
    ssf(i) = (yf(i+1)-yf(i))/Dsf(i);
    csf(i) = (xf(i+1)-xf(i))/Dsf(i);
end

Ax = zeros(ri-1,rf-1);
Ay = zeros(ri-1,rf-1);   
An = zeros(ri-1,rf-1); 
At = zeros(ri-1,rf-1); 
Bx = zeros(ri-1,rf-1);
By = zeros(ri-1,rf-1);   
Bn = zeros(ri-1,rf-1); 
Bt = zeros(ri-1,rf-1); 

for i=1:1:rf-1
    xp = [];
    yp = [];

    for j=1:1:ri
        xp2 = xi(j) - midx(i);
        yp2 = yi(j) - midy(i);
        xp(j,1) = [csf(i) ssf(i)]*[xp2;yp2];
        yp(j,1) = [-ssf(i) csf(i)]*[xp2;yp2];
    end

    for j=1:1:ri-1
        Vt = 0;
        Vn = 0;
        Vtg = 0;
        Vng = 0;
        
        midxp = (xp(j+1)+xp(j))/2;
        midyp = (yp(j+1)+yp(j))/2;
        nor = 4*pi;%Dsf(i)*2;
        num = (midxp + Dsf(i)/2)^2+midyp^2;
        den = (midxp - Dsf(i)/2)^2+midyp^2;
        Vt = log(num/den)/(nor);
        Vng = -Vt;
        
        num = midyp*Dsf(i);
        den = (midxp^2 + midyp^2 - (Dsf(i)/2)^2);
        Vn = 2*atan2(num,den)/(nor);
        Vtg = Vn;
        if ri == rf && xi(1) == xf(1) && xi(end) == xf(end) && i == j
            Vn = 0.5;
            Vt = 0;
            Vng = 0;
            Vtg = 0.5;
        end
            
        Ax(j,i) = Vt*csf(i) - Vn*ssf(i);
        Ay(j,i) = Vt*ssf(i) + Vn*csf(i);
        Bx(j,i) = Vtg*csf(i) - Vng*ssf(i);
        By(j,i) = Vtg*ssf(i) + Vng*csf(i);
        An(j,i) = -Ax(j,i)*ss(j) + Ay(j,i)*cs(j);
        At(j,i) = Ax(j,i)*cs(j) + Ay(j,i)*ss(j);
        Bn(j,i) = -Bx(j,i)*ss(j) + By(j,i)*cs(j);
        Bt(j,i) = Bx(j,i)*cs(j) + By(j,i)*ss(j); 

    end
end

if ri == rf && xi(1) == xf(1) && xi(end) == xf(end)
    da = diag(An);
    D1 = diag(da);
    An = An - D1;
    D2 = 0.5*ones(ri-1,1);
    D3 = diag(D2);
    An = An + D3;
    D1 = [];
    D2 = [];
    D3 = [];
    da = [];
    da = diag(At);
    D1 = diag(da);
    At = At - D1;
    da = [];
    D1 = [];
    da = diag(Bt);
    D1 = diag(da);
    Bt = Bt - D1;
    da = [];
    D1 = [];
    da = diag(Bn);
    D1 = diag(da);
    Bn = Bn - D1;
    Bt = Bt-diag(diag(Bt))+diag(0.5*ones(length(At),1));
    

end       
        
        
        
    


end

