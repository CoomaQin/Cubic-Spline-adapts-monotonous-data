%improved cubic spline 
%x are nodes, y are data. dy0,dyn are derivations of left and right
%endpoints. xx is interpolated interval.
function s=monocubicspline(x,y,dy0,dyn,xx)
n=length(x)-1;
h=diff(x);lamda=h(2:n)./(h(1:n-1)+h(2:n));mu=1-lamda;
g=3*(lamda.*diff(y(1:n))./h(1:n-1)+mu.*diff(y(2:n+1))./h(2:n));
delta=(y(2:n+1)-y(1:n))./h(1:n);
%base functions
function y=alpha0(x)
y=2*x.^3-3*x.^2+1;
end
function y=alpha1(x)
y=-2*x.^3+3*x.^2;
end
function y=beta0(x)
y=x.^3-2*x.^2+x;
end
function y=beta1(x)
y=x.^3-x.^2;
end
function y=phi(a,b)
    y=a-(2*a+b-3).^2./(3*(a+b-2));
end
%boundary conditions
g(1)=g(1)-lamda(1)*dy0;
g(n-1)=g(n-1)-mu(n-1)*dyn;
dy=nachase(lamda,2*ones(1:n-1),mu,g);
m=[dy0;dy;dyn];
for i=2:n-1
    if(m(i)<0)
        m(i)=delta(i-1)+h(i).*(delta(i)-delta(i-1))./(x(i+1)-x(i-1));
    end
end

if nargin>=5
    s=zeros(size(xx));
    for i=1:n
        if i==1
            kk=find(xx<=x(2));
        elseif i==n
            kk=find(xx>x(n));
        else
            kk=find(xx>x(i)&xx<=x(i+1));
        end
        xbar=(xx(kk)-x(i))/h(i);
        s(kk)=alpha0(xbar)*y(i)+alpha1(xbar)*y(i+1)+h(i)*beta0(xbar)*m(i)+h(i)*beta1(xbar)*m(i+1);
    end
end

end

