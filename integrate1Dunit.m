function integrate1Dunit(n)
%DESCRIPTIONS:
% Compute single int from a to b with n collocation points. -1<=a <= b<=1.

a = -1;
b = 1;
if(n<1)
    error('Need more points to evaluate integrals!');
end

%Integrate these test functions:
% f = @(x) x.^2;
f = @(x) 1./x;
% f = @(x) 1./x.^2;%does not conv
% f = @(x) 1./x.^3;%Grad Green kernel
% f = @(x) log(x);
% f = @(x) (x.^2)./sqrt(1-x.^2); %Chebyshev1
% f = @(x) 1./sqrt(1-x.^2);
% f = @(x) 1./(1-x.^2);%does not conv
% f = @(x) (x.^2).*exp(-x);

%Define the standard weight and points and change variable:
x = linspace(a,b,n);
w = (b-a)/(n-1);

%Call different algorithms to compare:
Matlab = integral(f,a,b)
% Collo = collocation(f)
% Trape = trapezoidal(f)
% Simps = Simpson(f)
QuadC1 = QuadChev1(f) %good for singular ints, but 1/x^2, 1/(1-x^2)
QuadC2 = QuadChev2(f) %good for singular ints even with only 10 collocation
%points but not stable (e.g. 1/x, 1/x^3 with 100 points), but 1/x^2, 1/(1-x^2)
%since these diverge.
ColloC1 = ColloChev1(f) %good for singular ints, but 1/x^2, 1/(1-x^2)
ColloC2 = ColloChev2(f) %comparable to ColloChev1

% GaussDR = Gaussian('dr') %nearly the best!!!
% GaussEK = Gaussian('ek')
% GaussSS = Gaussian('ss')
% GaussLO = Gaussian('lo')
% GaussRA = Gaussian('ra')
% GaussNCC = Gaussian('ncc')
% GaussNCO = Gaussian('nco')
% GaussNC = Gaussian('nc')

% f = @(x) x.^2; %for this original function: f = @(x) (x.^2)./sqrt(1-x.^2);
% GaussC = GaussianChev1()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function S = collocation(f)
        %The integral:
        %Integral (-1 <= X <= 1) F(X)dX
        %The quadrature rule:
        %Sum (1 <= I <= N) W(I)*F(X(I))
        
        y = f(x);
        S = sum(y*w);
    end

    function S = trapezoidal(f)
        %The integral:
        %Integral (-1 <= X <= 1) F(X)dX
        %The quadrature rule:
        %Sum (1 <= I <= N) W(I)*F(X(I))
        
        y = f(x);
        y(1) = y(1)/2;
        y(n) = y(n)/2;
        S = sum(y*w);
        
        %y = f(x);
        %S = sum(diff(x).*(y(1:end-1)+y(2:end))/2);
    end

    function S = Simpson(f)
        %The integral:
        %Integral (-1 <= X <= 1) F(X)dX
        %The quadrature rule:
        %Sum (1 <= I <= N) W(I)*F(X(I))
        
        y = 4*f(x);
        y(1:2:n) = y(1:2:n)/2;
        y(1) = f(x(1));
        y(n) = f(x(n));
        S = sum(y*w/3);
    end

    function S = Gaussian(alg)
        %The integral:
        %Integral (-1 <= X <= 1) F(X)dX
        %The quadrature rule:
        %Sum (1 <= I <= N) W(I)*F(X(I))
        
        xc = x;
        switch alg
            case 'dr'
                [xc,wc] = legendre_dr_compute(n); %compute roots of Legendre poly and weights
            case 'ek'
                [xc,wc] = legendre_ek_compute(n); %compute roots of Legendre poly and weights
            case 'ss'
                [xc,wc] = legendre_ss_compute(n); %compute roots of Legendre poly and weights
            case 'lo'
                [xc,wc] = lobatto_compute(n);
            case 'nc'
                wc = nc_compute_weights(n,x(1),x(n),x); %Newton-Cotes quadrature rule
                wc = wc.';
            case 'ncc'
                [xc,wc] = ncc_compute(n);%Newton-Cotes close
            case 'nco'
                [xc,wc] = nco_compute(n);%Newton-Cotes open
            otherwise
                [xc,wc] = radau_compute(n);
        end
        
        y = f(xc);
        S = sum(y.*wc);
    end

    function S = GaussianChev1()
        %The integral:
        %Integral (-1 <= X <= 1) F(X) / sqrt (1 - x^2) dX
        %The quadrature rule:
        %Sum (1 <= I <= N) W(I)*F(X(I))
        
        xc = cos((2*n+1-2*(1:n))*pi/(2*n)); %Chebyshev roots
        wc = pi/n; %Chebyshev weight
        
        y = f(xc);
        S = sum(y*wc);
    end

    function S = QuadChev1(f)
        %The integral:
        %Integral (-1 <= X <= 1) F(X)dX
        %The quadrature rule using Chebychev type 1:
        %Sum (1 <= I <= N) W(I)*F(X(I))
        
        A = zeros(n);
        RHS = zeros(n,1);
        %xc = cos(pi*(1+2*(0:n-1))/(2*(n+1))); %Chebyshev roots of Tn 1st
        %kind: incorrect result for 1/x.
        xc = cos(pi*(1:n)/(n+1)); %Chebyshev roots of Un 2nd kind
        for i=1:n
            fun = @(x) cos(i.*acos(x)); %Chebyshev Tn 1st kind
            A(i,:) = fun(xc);
            RHS(i) = integral(fun,a,b);
        end
        wc = gmres(A,RHS);
        
        y = f(xc);
        S = sum(y*wc);
    end

    function S = QuadChev2(f)
        %The integral:
        %Integral (-1 <= X <= 1) F(X)dX
        %The quadrature rule using Chebychev type 2:
        %Sum (1 <= I <= N) W(I)*F(X(I))
        
        A = zeros(n);
        RHS = zeros(n,1);
        xc = cos(pi*(1:n)/(n+1)); %Chebyshev roots of Un 2nd kind
        for i=1:n
            fun = @(x) sin(i.*acos(x))./sin(acos(x)); %Chebyshev Un 2nd kind
            A(i,:) = fun(xc);
            RHS(i) = integral(fun,a,b);
        end
        wc = gmres(A,RHS);
        
        y = f(xc);
        S = sum(y*wc);
    end

    function S = ColloChev1(f)
        %Computing an integral using Chebychev approximation type 1:
        %Integral (-1 <= X <= 1) F(X)dX
        
        %xc = cos(pi*(1+2*(0:n))/(2*(n+1))); %Chebyshev roots of Tn 1st
        %kind: incorrect result.
        xc = cos(pi*(1:n)/(n+1)); %Chebyshev roots of Un 2nd kind
        c = zeros(1,n);
        F = f(xc);
        AC = acos(xc);
        for k=0:n
            c(k+1) = sum(F.*cos(k*AC));
        end
        c = c/(n+1);
        c(2:n+1) = 2*c(2:n+1);
        S = 0;
        for k =0:n
            Tk = @(x) cos(k.*acos(x));%Chebyshev function Tn 1st kind
            S = S+c(k+1)*integral(Tk,a,b);
        end
    end

    function S = ColloChev2(f)
        %Computing an integral using Chebychev approximation type 2:
        %Integral (-1 <= X <= 1) F(X)dX
        
        %xc = cos(pi*(1+2*(0:n))/(2*(n+1))); %Chebyshev roots of Tn 1st kind: incorrect
        xc = cos(pi*(1:n)/(n+1)); %Chebyshev roots of Un 2nd kind
        c = zeros(1,n);
        F = f(xc);
        AC = acos(xc);
        for k=1:n
            c(k) = sum((1-xc.^2).*F.*sin(k*AC)./sin(AC));
        end
        c = 2/n*c;
        S = 0;
        for k =1:n
            Uk = @(x) sin(k.*acos(x))./sin(acos(x)); %Chebyshev Un 2nd kind
            S = S+c(k)*integral(Uk,a,b);
        end
    end

end