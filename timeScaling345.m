function [s, sdot, sddot] = timeScaling345(l_all, T, t)

Phi = t(:)'/T;

s = l_all*(6*Phi.^5 - 15*Phi.^4 + 10*Phi.^3);

sdot  = (l_all/T)  * (30*Phi.^4 - 60*Phi.^3 + 30*Phi.^2);
sddot = (l_all/T^2)* (120*Phi.^3 - 180*Phi.^2 + 60*Phi);

end
