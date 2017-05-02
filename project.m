clear
clc

%% Initiallization 
N=10; % Number of nodes in x or y direction
h=2*pi/(N-1); % step size
x=linspace(-pi-h,pi,N+1); 
y=linspace(pi,-pi,N);
u_domain=zeros(N,N+1); % domain has additional nodes in x to use later as ghost nodes
[X,Y]=meshgrid(x,y); % x,y grid for plotting

%% Domain Boundary Conditions
% Top B.C.'s
f=x(2:end).*(x(2:end)+pi).^2;
u_domain(1,2:end)=f;

% Bottom B.C.'s
g=(x(2:end)+pi).^2.*cos(-x(2:end));
u_domain(N,2:end)=g;

% Right B.C.'s
r=g(end)+((y+pi)/(2*pi))*(f(end)-g(end));
u_domain(:,N+1)=r';

clear r g f

%% Forcing Function 
for k=1:length(y)
    for j=1:length(x)
        Force(k,j)=(sin(pi*(x(j)+pi)/(2*pi))*cos((pi/2)*(2*(y(k)+pi)/(2*pi)+1)));
    end
end

%% Gauss Seidel
u=u_domain;
gauss_errorval=100;
gauss_iterations=0;
while gauss_errorval>1
    gauss_iterations=gauss_iterations+1;
    
    u_old=u; % This value is used to calculate error from iteration to iteration
    
    % The while loop below solves the values in the domain by starting at
    % the outermost unknowns near the direchlet boundary conditions and
    % working towards the center of the domain, utilizing as many solved
    % values as possible to reduce the necessary number of iterations.
    counter=2;
    while counter<=N/2
        
        % Top Pyramid
        for j=      N-counter+2   :  -1 :   counter+1
            u( counter , j ) = (Force(counter,j)*h^2+...
                u(counter-1,j)+u(counter+1,j)+u(counter,j-1)+u(counter,j+1))*0.25;
        end
        
        % Bottom Pyramid
        for j=      N-counter+2  :  -1 :    counter+1
            u( N-counter+1 , j ) = (Force(N-counter+1,j) *h^2+...
                u(N-counter+2,j)+u(N-counter,j)+u(N-counter+1,j-1)+u(N-counter+1,j+1))*0.25;
        end
        
        % Right Pyramid
        for k=       counter+1      :        N-counter
            u( k , N-counter+2 ) = (Force(k,N-counter+2)*h^2+...
                u(k-1,N-counter+2)+u(k+1,N-counter+2)+u(k,N-counter+1)+u(k,N-counter+3))*0.25;
        end
        
        counter=counter+1;
    end
    
    % Left Pyramid
    counter2=0;
    if mod(N,2)==0 % even number of nodes
        for j=       floor(N/2)+1 :   -1 :  3
            counter2=counter2+1;
            for k =  floor(N/2)-counter2+2    :    floor(N/2)+counter2-1
                u( k , j) = (Force(k,j)*h^2+...
                    u(k-1,j)+u(k+1,j)+u(k,j-1)+u(k,j+1))*0.25;
            end
            
        end
        
    else % odd number of nodes
        for j=       floor(N/2)+2 :   -1 :  3
            counter2=counter2+1;
            for k =  floor(N/2)-counter2+2   :    floor(N/2)+counter2
                u( k , j) = (Force(k,j)*h^2+...
                    u(k-1,j)+u(k+1,j)+u(k,j-1)+u(k,j+1))*0.25;
            end
            
        end
    end
    
    % ghost nodes; this code duplicates values from the right side of the
    % Nuemann boundary to the left side, in order to simulate the du/dx=0
    % condition
    u(:,1)=u(:,3);
    
    % left boundary nodes
    for k=2:N-1
        u(k,2)=(Force(k,2)*h^2+...
            u(k-1,2)+u(k+1,2)+u(k,1)+u(k,3))*0.25;
    end
    
    % remaining 2 nodes at nueman boundary
    u(1,2)=(Force(1,2)*h^2+...
        u(2,2)+u(1,1)+u(1,3))*0.25;
    u(N,2)=(Force(N,2)*h^2+...
        u(N-1,2)+u(N,1)+u(N,3))*0.25;
    
    for k=1:N
        for j=1:N+1
            if u(k,j) ~= 0
                error(k,j)=abs((u(k,j)-u_old(k,j))/u(k,j))*100;
            end
        end
    end
    % average error over domain
    gauss_errorval=mean(mean(error));
end
gauss_iterations
gauss_errorval
figure(1)
mesh(X,Y,u),xlabel('x'),ylabel('y'),zlabel('u'),title('Gauss Seidel Method')
clear error

%% Gauss-Seidel with Relaxation

u_relaxed=u_domain;

errorval_relaxed=100;
iterations_relaxed=0;
lambda=1.4;
while errorval_relaxed>1
    iterations_relaxed=iterations_relaxed+1;
    u_old2=u_relaxed; % This value is used to calculate error from iteration to iteration
    % Circling
    counter=2;
    while counter<=N/2
        
        % Top Pyramid
        for j=      N-counter+2   :  -1 :   counter+1
            u_relaxed( counter , j ) = ((sin(pi*(x(j)+pi)/(2*pi))*cos((pi/2)*(2*(y(counter)+pi)/(2*pi)+1)))*h^2+...
                u_relaxed(counter-1,j)+u_relaxed(counter+1,j)+u_relaxed(counter,j-1)+u_relaxed(counter,j+1))*0.25;
        end
        
        
        % Bottom Pyramid
        for j=      N-counter+2  :  -1 :    counter+1
            u_relaxed( N-counter+1 , j ) = ((sin(pi*(x(j)+pi)/(2*pi))*cos((pi/2)*(2*(y(N-counter+1)+pi)/(2*pi)+1)))*h^2+...
                u_relaxed(N-counter+2,j)+u_relaxed(N-counter,j)+u_relaxed(N-counter+1,j-1)+u_relaxed(N-counter+1,j+1))*0.25;
        end
        
        % Right Pyramid
        for k=       counter+1      :        N-counter
            u_relaxed( k , N-counter+2 ) = ((sin(pi*(x(N-counter+2)+pi)/(2*pi))*cos((pi/2)*(2*(y(k)+pi)/(2*pi)+1)))*h^2+...
                u_relaxed(k-1,N-counter+2)+u_relaxed(k+1,N-counter+2)+u_relaxed(k,N-counter+1)+u_relaxed(k,N-counter+3))*0.25;
        end
        
        
        counter=counter+1;
        
    end
    
    % Left Pyramid
    counter2=0;
    if mod(N,2)==0 % even number of nodes
        for j=       floor(N/2)+1 :   -1 :  3
            counter2=counter2+1;
            for k =  floor(N/2)-counter2+2    :    floor(N/2)+counter2-1
                u_relaxed( k , j) = ((sin(pi*(x(j)+pi)/(2*pi))*cos((pi/2)*(2*(y(k)+pi)/(2*pi)+1)))*h^2+...
                    u_relaxed(k-1,j)+u_relaxed(k+1,j)+u_relaxed(k,j-1)+u_relaxed(k,j+1))*0.25;
            end
            
        end
        
    else % odd number of nodes
        for j=       floor(N/2)+2 :   -1 :  3
            counter2=counter2+1;
            for k =  floor(N/2)-counter2+2   :    floor(N/2)+counter2
                u_relaxed( k , j) = ((sin(pi*(x(j)+pi)/(2*pi))*cos((pi/2)*(2*(y(k)+pi)/(2*pi)+1)))*h^2+...
                    u_relaxed(k-1,j)+u_relaxed(k+1,j)+u_relaxed(k,j-1)+u_relaxed(k,j+1))*0.25;
            end
            
        end
    end
    
    % ghost nodes
    u_relaxed(:,1)=u_relaxed(:,3);
    
    % left boundary nodes
    for k=2:N-1
        u_relaxed(k,2)=((sin(pi*(x(2)+pi)/(2*pi))*cos((pi/2)*(2*(y(k)+pi)/(2*pi)+1)))*h^2+...
            u_relaxed(k-1,2)+u_relaxed(k+1,2)+u_relaxed(k,1)+u_relaxed(k,3))*0.25;
    end
    
    % remaining 2 nodes at nueman boundary
    u_relaxed(1,2)=((sin(pi*(x(2)+pi)/(2*pi))*cos((pi/2)*(2*(y(1)+pi)/(2*pi)+1)))*h^2+...
        u_relaxed(2,2)+u_relaxed(1,1)+u_relaxed(1,3))*0.25;
    u_relaxed(N,2)=((sin(pi*(x(2)+pi)/(2*pi))*cos((pi/2)*(2*(y(N)+pi)/(2*pi)+1)))*h^2+...
        u_relaxed(N-1,2)+u_relaxed(N,1)+u_relaxed(N,3))*0.25;
    
    
    for k=1:N
        for j=1:N+1
            error(k,j)=abs((u_relaxed(k,j)-u_old2(k,j))/u_relaxed(k,j))*100;
        end
    end
    % average error over domain
    errorval_relaxed=mean(mean(error));
    u_relaxed=lambda*u_relaxed+(1-lambda)*u_old2;
end
iterations_relaxed
errorval_relaxed
figure(2)
mesh(X,Y,u_relaxed),xlabel('x'),ylabel('y'),zlabel('u'),title('Gauss Seidel with Relaxation')
clear error
%% Gauss Siedel No Forcing Function

% Boundary Conditions
uNoF=u_domain;

errorval_no_force=100;
iterations_no_force=0;
while errorval_no_force>1
    iterations_no_force=iterations_no_force+1;
    u_oldNoF=uNoF; % This value is used to calculate error from iteration to iteration
    
    
    counter=2;
    while counter<=N/2
        
        % Top Pyramid
        for j=      N-counter+2   :  -1 :   counter+1
            uNoF( counter , j ) =   (uNoF(counter-1,j)+uNoF(counter+1,j)+uNoF(counter,j-1)+uNoF(counter,j+1))*0.25;
        end
        
        
        % Bottom Pyramid
        for j=      N-counter+2  :  -1 :    counter+1
            uNoF( N-counter+1 , j ) =   (uNoF(N-counter+2,j)+uNoF(N-counter,j)+uNoF(N-counter+1,j-1)+uNoF(N-counter+1,j+1))*0.25;
        end
        
        % Right Pyramid
        for k=       counter+1      :        N-counter
            uNoF( k , N-counter+2 ) =   (uNoF(k-1,N-counter+2)+uNoF(k+1,N-counter+2)+uNoF(k,N-counter+1)+uNoF(k,N-counter+3))*0.25;
        end
        
        
        counter=counter+1;
        
    end
    
    % Left Pyramid
    counter2=0;
    if mod(N,2)==0 % even number of nodes
        for j=       floor(N/2)+1 :   -1 :  3
            counter2=counter2+1;
            for k =  floor(N/2)-counter2+2    :    floor(N/2)+counter2-1
                uNoF( k , j) =        (uNoF(k-1,j)+uNoF(k+1,j)+uNoF(k,j-1)+uNoF(k,j+1))*0.25;
            end
            
        end
        
    else % odd number of nodes
        for j=       floor(N/2)+2 :   -1 :  3
            counter2=counter2+1;
            for k =  floor(N/2)-counter2+2   :    floor(N/2)+counter2
                uNoF( k , j) =        ( uNoF(k-1,j)+uNoF(k+1,j)+uNoF(k,j-1)+uNoF(k,j+1))*0.25;
            end
            
        end
    end
    
    % ghost nodes  
    uNoF(:,1)=uNoF(:,3);
    
    % left boundary nodes
    for k=2:N-1
        uNoF(k,2)=           (uNoF(k-1,2)+uNoF(k+1,2)+uNoF(k,1)+uNoF(k,3))*0.25;
    end
    
    % remaining 2 nodes at nueman boundary
    uNoF(1,2)=         ( uNoF(2,2)+uNoF(1,1)+uNoF(1,3))*0.25;
    uNoF(N,2)=         ( uNoF(N-1,2)+uNoF(N,1)+uNoF(N,3))*0.25;
    
    
    for k=1:N
        for j=1:N+1
            error3(k,j)=abs((uNoF(k,j)-u_oldNoF(k,j))/uNoF(k,j))*100;
        end
    end
    % average error over domain
    errorval_no_force=mean(mean(error3));
end
iterations_no_force
errorval_no_force
figure(3)
mesh(X,Y,uNoF),xlabel('x'),ylabel('y'),zlabel('u'),title('Gauss Seidel with no Forcing Function')

clear u_old u_old2 u_oldNoF u_domain counter counter2 j k error_no_force x y
