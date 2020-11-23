%% Question 5

%Alex Vanderhoeff
%250858436
%For Phys3330

%relax - program to solve laplace equation 
clear all
%initialize parameters
method = menu('Numerical method','Jacobi');
N=input('enter number of grid points on a side: ');
L=100; %system size
h=L/(N-1); %grid spacing
x=(0:N-1)*h; %x coord
y=(0:N-1)*h; %y coord

% set initial guess as first term in separation of variables soln
phi0=1; %pot at y=1
phi=phi0*4/(pi*sinh(pi))*sin(pi*x'/L)*sinh(pi*y/L);

% set Boundary conditions
phi(1,:)=y/100;
phi(N,:)=0;
phi(:,1)=sin((pi*x)/100);
phi(:,N)=cos((3*pi*x)/100);
% fprintf('Potential at y=L equals %g \n', phi0);
% fprintf('Potential is zero

%Loop until number iterations wanted achieved
flops(0); %reset number of flops to 0
newphi=phi; %copy solution
itermax=N^2; %avoids overly long runs
changedesired=1e-4; %stop when change is given fraction
fprintf('Desired fractional change = %g\n',changedesired);
for iter=1:itermax
    changesum=0;
    if(method==1) %Jacobi method
        for i=2:(N-1) %loop over interior points
            for j=2:(N-1)
                newphi(i,j)=.25*(phi(i+1,j)+phi(i-1,j)+phi(i,j-1)+phi(i,j+1));
                changesum=changesum+abs(1-phi(i,j)/newphi(i,j));
            end
        end
        phi=newphi;
    end


%check if fractional change small enough to stop
    change(iter)=changesum/(N-2)^2;
    if(rem(iter,10) < 1)
        fprintf('After %g iterations, fractional change = %g\n', iter,change(iter));
    end
    if(change(iter) < changedesired)
        fprintf('Desired accuracy achieved after %g iterations\n',iter);
        fprintf('Breaking out of main loop\n');
        break;
    end
end

%plot potential as contour and surface plots

figure(1);
clf;
clevels=0:(0.1):1; %contour levels
cs=contour(x,y,flipud(rot90(phi)),clevels);
xlabel('X');
ylabel('Y');
clabel(cs);
title(sprintf('Potential after %g iterations',iter));

figure(2);
clf;
mesh(x,y,flipud(rot90(phi)));
xlabel('X');
ylabel('Y');
zlabel('\phi(x,y)');

%plot fractional change vs iteration
figure(3);
clf;
semilogy(change);
xlabel('Iteration');
ylabel('Fractional change');
title(sprintf('Number of flops = %g\n',flops));




