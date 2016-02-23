fun = @(x,y) (-1/4*pi)*sin(y);
q0 = integral2(fun,0,pi,0,2*pi)


fun = @(x,y) 1./( sqrt(x + y) .* (1 + x + y).^2 );
polarfun = @(theta,r) fun(r.*cos(theta),r.*sin(theta)).*r;
rmax = @(theta) 1./(sin(theta) + cos(theta));
q1 = integral2(polarfun,0,pi/2,0,rmax)


fun = @(x,y,z) -1./((x.^2+y.^2+z.^2).^3).*(x.^2)*(1/(4*pi));
r = 1;
polarfun = @(theta,phi) fun(r.*cos(theta).*sin(phi),r.*sin(theta).*sin(phi),r.*cos(phi));
q2 = integral2(polarfun,0,2*pi,0,pi)


fun = @(s,t) -(s-t)./((abs(s-t).^3)*(4*pi));
r = 1;
polarfun = @(theta1,phi1,theta2,phi2) fun([r.*cos(theta1).*sin(phi1),r.*sin(theta1).*sin(phi1),r.*cos(phi1)],[r.*cos(theta2).*sin(phi2),r.*sin(theta2).*sin(phi2),r.*cos(phi2)]);
q3 =  integral2(integral2(polarfun,0,2*pi,0,pi,theta2,phi2),0,2*pi,0,pi)
