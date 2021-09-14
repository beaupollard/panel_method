function [phi,cp,Pressure] = Determine_Forces(Vt,Vn,n,t,pos,m,gamma,midx...
    ,midy,vortx,vorty,del,x_wp,y_wp,delta,phi,Vdx,Vdy,x,y,gamma_w,Vbt,Vb)
%Determine the forces acting on the body

global k
global interval

theta = pos(k,3);
Vx = Vt.*t(:,1)+Vn.*n(:,1);
Vy = Vt.*t(:,2)+Vn.*n(:,2);

%% Velocities with respect to body fixed frame
ux = Vx.*cos(theta) - Vy.*sin(theta);
uy = Vx.*sin(theta) + Vy.*cos(theta);

if k > 100
    avg_theta = mean(pos(k-100:k,3));
else
    avg_theta = 0;
end

%% Compute Phi
[phinew] = Potential2(m,gamma,midx,midy,vortx,vorty,del,n,x,y...
        ,gamma_w,x_wp,y_wp,delta,Vt,avg_theta,k,pos(k,:),Vb);

phi(:,k) = phinew;
dphi = (phi(:,k)-phi(:,k-1))./interval;

for i=1:length(Vx)
    cp(i,1) = (-1/2*(ux(i)^2+uy(i)^2) - (...
        -Vdx(i)*ux(i)-Vdy(i)*uy(i)) - dphi(i));
end
Pressure = cp;


end

