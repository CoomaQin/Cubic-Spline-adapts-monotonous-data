function m=AKIMA(x,y,dy0,dyn,xx)
%   AKIMA algorithm
%x are nodes, y are data. dy0,dyn are derivations of left and right
%endpoints. xx is interpolated interval.
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
%boundary conditions
g(1)=g(1)-lamda(1)*dy0;
g(n-1)=g(n-1)-mu(n-1)*dyn;
dy=nachase(lamda,2*ones(1:n-1),mu,g);
m=[dy0;dy;dyn];
a(3:n-2)=abs(delta(4:n-1)-delta(3:n-2));
b(3:n-2)=abs(delta(2:n-3)-delta(1:n-4));
m(3:n-2)=a(3:n-2).*delta(2:n-3)./(a(3:n-2)+b(3:n-2))+b(3:n-2).*delta(3:n-2)./(a(3:n-2)+b(3:n-2));
for i=1:n
    if(m(i)<0)
        m(i)=-m(i);
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