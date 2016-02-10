function integrate1D(n,a,b)
%DESCRIPTIONS: 
%   Compute single int from a to b with n collocation points.

if(n<1)
    error('Need more points to evaluate integrals!');
end

%Integrate these test functions:
%       f = @(x) x.^2;
%       f = @(x) 1./x;
%       f = @(x) 1./x.^2;              %does not conv
       f = @(x) 1./x.^3;              %Grad Green kernel
%       f = @(x) log(x);
%       f = @(x) (x.^2)./sqrt(1-x.^2); %Chebyshev1
%       f = @(x) 1./sqrt(1-x.^2);
%       f = @(x) 1./(1-x.^2);          %does not conv
%       f = @(x) (x.^2).*exp(-x);

%Define the standard weight and points and change variable:
    x = linspace(a,b,n);
    w = (b-a)/(n-1);    

%Call different algorithms to compare:    
    Matlab = integral(f,a,b)
%     Collo = collocation(f)
%     Trape = trapezoidal(f)
%     Simps = Simpson(f)
    QuadC1 = QuadChev1(f)   %good for singular ints, but 1/x^2, 1/(1-x^2), 1/x, 1/x^3
    QuadC2 = QuadChev2(f)   %conv better than QuadChev1(), works for 1/x, 1/x^3 (10 points) 
% These two are much reliable and better than QuadChev1(), QuadChev2():     
    ColloC1 = ColloChev1(f) %good for singular ints, but 1/x^2, 1/(1-x^2): diverge
    ColloC2 = ColloChev2(f) %comparable to ColloChev1()
    
    LG = LegendreGauss(f)
    
    GaussDR = Gaussian('dr')   %nearly the best!!!
%     GaussEK = Gaussian('ek') 
%     GaussSS = Gaussian('ss') 
%     GaussLO = Gaussian('lo') 
%     GaussRA = Gaussian('ra')
%     GaussNCC = Gaussian('ncc') 
%     GaussNCO = Gaussian('nco')
%     GaussNC = Gaussian('nc') 
    
%     f = @(x) x.^2; %for this original function: f = @(x) (x.^2)./sqrt(1-x.^2);
%     GaussC = GaussianChev1()    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function S = collocation(f)
%    The integral:
%      Integral (-1 <= X <= 1) F(X)dX
%    The quadrature rule:
%      Sum (1 <= I <= N) W(I)*F(X(I))

        y = f(x);        
        S = sum(y*w);
    end

    function S = trapezoidal(f)
%    The integral:
%      Integral (-1 <= X <= 1) F(X)dX
%    The quadrature rule:
%      Sum (1 <= I <= N) W(I)*F(X(I))

        y = f(x);
        y(1) = y(1)/2;
        y(n) = y(n)/2;        
        S = sum(y*w);
        
        %y = f(x);
        %S = sum(diff(x).*(y(1:end-1)+y(2:end))/2);
    end

    function S = Simpson(f)
%    The integral:
%      Integral (-1 <= X <= 1) F(X)dX
%    The quadrature rule:
%      Sum (1 <= I <= N) W(I)*F(X(I))

        y = 4*f(x);
        y(1:2:n) = y(1:2:n)/2;
        y(1) = f(x(1));
        y(n) = f(x(n));        
        S = sum(y*w/3);
    end

    function S = Gaussian(alg) 
%    The integral:
%      Integral (-1 <= X <= 1) F(X)dX
%    The quadrature rule:
%      Sum (1 <= I <= N) W(I)*F(X(I))

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
                [xc,wc] = ncc_compute(n);  %Newton-Cotes close                
            case 'nco'
                [xc,wc] = nco_compute(n);  %Newton-Cotes open             
            otherwise
                [xc,wc] = radau_compute(n);               
        end    
        
        y = f(xc);  
        S = sum(y.*wc);
    end

    function S = GaussianChev1() 
%    The integral:
%      Integral (-1 <= X <= 1) F(X) / sqrt (1 - x^2) dX
%    The quadrature rule:
%      Sum (1 <= I <= N) W(I)*F(X(I))

        xc = cos((2*n+1-2*(1:n))*pi/(2*n)); %Chebyshev roots   
        wc = pi/n; %Chebyshev weight
        
        y = f(xc);       
        S = sum(y*wc);
    end

    function S = QuadChev1(f) 
%    The integral:
%      Integral (-1 <= X <= 1) F(X)dX
%    The quadrature rule using Chebychev type 1:
%      Sum (1 <= I <= N) W(I)*F(X(I))

        A = zeros(n);
        RHS = zeros(n,1);
        J = (pi*(1:n)/(n+1));
        yc = cos(J); %Chebyshev roots of U(n+1), 2nd kind, or extrema points of T(n+1),
        %1st kind, (without the first point when k=0)
        xc = yc*(b-a)/2+(b+a)/2; %Change of variable: yc in [-1,1] -> xc in [a,b]        
        %xc = cos(pi*(1+2*(0:n-1))/(2*(n+1))); %Chebyshev roots of Tn 1st
        %kind: incorrect result for 1/x.
        for k=0:2:n-1 %must start from 0
            fun = @(t) cos(k.*t); %Chebyshev Tn 1st kind
            A(k+1,:) = fun(J);
            %fun = @(x) cos(k.*acos(x)); %Chebyshev Tn 1st kind
            %RHS(k+1) = integral(fun,-1,1);
            if(mod(k,2)==0)
                RHS(k+1) = 2/(1-k^2);
            end
        end          
        wc = gmres(A,RHS*(b-a)/2);          

        y = f(xc);            
        S = sum(y*wc);
    end

    function S = QuadChev2(f)
        %    The integral:
        %      Integral (-1 <= X <= 1) F(X)dX
        %    The quadrature rule using Chebychev type 2:
        %      Sum (1 <= I <= N) W(I)*F(X(I))
        
        A = zeros(n);
        RHS = zeros(n,1);
        J = (pi*(1:n)/(n+1));
        yc = cos(J); %Chebyshev roots of U(n+1), 2nd kind
        xc = yc*(b-a)/2+(b+a)/2; %Change of variable: yc in [-1,1] -> xc in [a,b]        
        for i=1:n %must start from 1
            fun = @(t) sin(i.*t)./sin(t); %Chebyshev Un 2nd kind
            A(i,:) = fun(J);
%             fun = @(x) sin(i.*acos(x))./sin(acos(x)); %Chebyshev Un 2nd kind            
%             RHS(i) = integral(fun,-1,1);
            RHS(i) = 2*(sin(pi*i/2)^2)/i;
        end        
        wc = gmres(A,RHS*(b-a)/2);       
        
        y = f(xc);
        S = sum(y*wc);
    end

    function S = ColloChev1(f)
        %    Computing an integral using Chebychev approximation T(n) type 1:
        %      Integral (-1 <= X <= 1) F(X)dX
        
        %xc = cos(pi*(1+2*(0:n))/(2*(n+1))); %Chebyshev roots of T(n+1) 1st
        %kind: incorrect result.        
        c = zeros(1,n+1);                
        J = (pi*(1:n)/(n+1));
        yc = cos(J); %Chebyshev roots of U(n+1), 2nd kind, or extrema points of T(n+1),
        %1st kind, (without the first point when k=0)
        xc = yc*(b-a)/2+(b+a)/2; %Change of variable: yc in [-1,1] -> xc in [a,b]
        F = f(xc);
        %Find the coefficients in Chebyshev approximation of f:
        for k=0:2:n %must start from 0         
            c(k+1) = sum(F.*cos(k*J));            
        end    
        c = c/n;
        c(2:n) = 2*c(2:n);
        S = 0;
        for k =0:2:n %must start from 0
%             Tk = @(x) cos(k.*acos(x)); %Chebyshev function T(k) 1st kind
%             S = S+c(k+1)*integral(Tk,-1,1);
            S = S+c(k+1)*2/(1-k^2);
        end
        S = S*(b-a)/2; %due to change of variable
    end

    function S = ColloChev2(f)
        %    Computing an integral using Chebychev approximation U(n) type 2:
        %      Integral (-1 <= X <= 1) F(X)dX
                        
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
%             Uk = @(x) sin(k.*acos(x))./sin(acos(x)); %Chebyshev Un 2nd kind
%             S = S+c(k)*integral(Uk,-1,1);
            S = S+c(k)*2*(sin(pi*k/2)^2)/k;
        end
        S = S*(b-a)/2; %due to change of variable
    end

    function S = LegendreGauss(f)
       [xc,wc] = lgwt(n,a,b); 
       y = f(xc);
       S = sum(y.*wc);
    end

end