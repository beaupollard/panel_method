function [Cx,Cy,Ct,Cn] = Free_Vort_Affect(vortx,vorty,points,gamma_w,n1,t1,pan)
%Find the affects of a point vortex with CW as positive circulation

global k
global panels

[r,c] = size(points);
r2 = length(vortx);
xp = points(:,1);
yp = points(:,2);

for i=1:r
    Vx = 0;
    Vy = 0;

    for j=1:k-2
        delx = xp(i)-vortx(j);
        dely = yp(i)-vorty(j);
        del = delx^2+dely^2;
        
        if r == length(vortx) && vortx(1) == xp(1) && i == j
            Vx = Vx;
            Vy = Vy;
        else
            Vx = Vx + gamma_w(j+1)*dely/(2*pi*del);
            Vy = Vy - gamma_w(j+1)*delx/(2*pi*del);
        end
    end
    Cx(i,1) = Vx;
    Cy(i,1) = Vy;
    
    if pan == 0;
        Cn(i,1) = [n1(i,1) n1(i,2)]*[Cx(i);Cy(i)];
        Ct(i,1) = [t1(i,1) t1(i,2)]*[Cx(i);Cy(i)];
    else
        Cn = 0;
        Ct = 0;
    end
end
    
ct2 = Ct;        
        



end

