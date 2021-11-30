function dGdX = dGdx(Udx,dt)


dGdX = - (Udx*dt/2 + eye(size(Udx,1)));