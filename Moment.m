function [M_b] = Moment(cp,del,n1,midx,midy,pos)
%Calculate moment on tail and body

Fx = -cp.*del.*n1(:,1)./2;
Fy = -cp.*del.*n1(:,2)./2;

M_b = 0;
for i=1:length(midx)
    M_b = M_b-(midy(i)-pos(1,2))*Fx(i) + (midx(i)-pos(1,1))*Fy(i);
end

end

