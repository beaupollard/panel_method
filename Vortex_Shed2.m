function [vort_posx,vort_posy] = Vortex_Shed2(vortx,vorty,gamma_w ...
    ,gamma,m,x_wp,y_wp,Vxw,Vyw,x,y,delta,theta,n1,alpha)
%Shed vorticity

global k
global V_inf
global interval
global L
n2 = 0;
t2 = 0;

xp2 = (x(end)+x_wp)/2;
yp2 = (y(end)+y_wp)/2;

if k == 2
    vort_posx(1,1) = xp2 + Vxw*interval;
    vort_posy(1,1) = yp2 + Vyw*interval;
else
    
    [Ax,Ay,Bx2,By2] = Panel_Vort(vortx,vorty,x,y);
    l2 = length(Bx2);
    for i=1:k-2
        Bx(i,1) = sum(Bx2(i,:));
        By(i,1) = sum(By2(i,:));
    end
    
    xf = [x(end);x_wp];
    yf = [y(end);y_wp];
    
    [XX,XX,Bwx,Bwy] = Panel_Vort(vortx,vorty,xf,yf);
    
    pan = 1;
    [Cx,Cy,XX,XX] = Free_Vort_Affect(vortx,vorty,[vortx vorty],gamma_w,n2,t2,pan);
    
    Vxp = Ax*m + gamma(k)*Bx + L/delta(k)*(gamma(k-1)-gamma(k)).*Bwx + Cx + V_inf(1);
    Vyp = Ay*m + gamma(k)*By + L/delta(k)*(gamma(k-1)-gamma(k)).*Bwy + Cy + V_inf(2);
    
    r = length(vortx);
    vort_posx(1:r,1) = vortx + Vxp*interval;
    vort_posy(1:r,1) = vorty + Vyp*interval;
    
   
    vort_posx(r+1,1) = xp2 + Vxw*interval;
    vort_posy(r+1,1) = yp2 + Vyw*interval;
end

end

