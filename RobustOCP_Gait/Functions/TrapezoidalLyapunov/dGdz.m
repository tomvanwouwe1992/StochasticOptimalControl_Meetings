function dGdZ = dGdz(Udz,dt)


dGdZ = eye(size(Udz,1)) - Udz*dt/2;