function integrate2D(n,a,b)

fun = @(x,y) 1./abs(x-y);
q0 = integral2(fun,0,pi,0,pi)

fun = @(x,y) 1./( sqrt(x + y) .* (1 + x + y).^2 );
polarfun = @(theta,r) fun(r.*cos(theta),r.*sin(theta)).*r;
rmax = @(theta) 1./(sin(theta) + cos(theta));
q1 = integral2(polarfun,0,pi/2,0,rmax)

fun = @(x,y,z) -1./((x.^2+y.^2+z.^2).^3).*(x.^2)*(1/(4*pi));
r = 1;
polarfun = @(theta,phi) fun(r.*cos(theta).*sin(phi),r.*sin(theta).*sin(phi),r.*cos(phi));
q2 = integral2(polarfun,0,2*pi,0,pi)

% fun = @(s,t) -(s-t)./((abs(s-t).^3)*(4*pi));
% r = 1;
% polarfun = @(theta1,phi1,theta2,phi2) fun([r.*cos(theta1).*sin(phi1),r.*sin(theta1).*sin(phi1),r.*cos(phi1)],[r.*cos(theta2).*sin(phi2),r.*sin(theta2).*sin(phi2),r.*cos(phi2)]);
% q3 =  integral2(integral2(polarfun,0,2*pi,0,pi,@theta2,@phi2),0,2*pi,0,pi)


    function S = ColloChev2(f)
        %Computing an integral using Chebychev approximation U(n) type 2:
        %  Integral (-1 <= X <= 1) F(X)dX
        
        c = zeros(1,n);
        J = (pi*(1:n)/(n+1));
        yc = cos(J); %Chebyshev roots of U(n+1), 2nd kind
        xc = yc*(b-a)/2+(b+a)/2; %Change of variable: yc in [-1,1] -> xc in [a,b]
        F = f(xc);
        %Find the coefficients in Chebyshev approximation of f:
        for k=1:2:n %must start from 1
            c(k) = sum((1-yc.^2).*F.*sin(k*J)./sin(J));
        end
        c = c*2/n;
        S = 0;
        for k =1:2:n %must start from 1
            % Uk = @(x) sin(k.*acos(x))./sin(acos(x)); %Chebyshev Un 2nd kind
            % S = S+c(k)*integral(Uk,-1,1);
            S = S+c(k)*2*(sin(pi*k/2)^2)/k;
        end
        S = S*(b-a)/2; %due to change of variable
    end

end