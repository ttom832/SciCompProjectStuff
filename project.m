
function[h,X,Y,cycles,percent_error,u_gauss_seidel,u_gauss_seidel_relaxed,u_gauss_seidel_no_force,u_gauss_seidel_no_force_relaxed]=project(Nodes,Lambda)


%% Initiallization
N=Nodes; % Number of nodes in x or y direction
h=2*pi/(N-1); % step size
x=linspace(-pi-h,pi,N+1); % x axis domain
y=linspace(pi,-pi,N); % y axis domain
u_domain=zeros(N,N+1); % domain has additional nodes in x to use later as ghost nodes
[X,Y]=meshgrid(x,y); % x,y grid for plotting
addon=Lambda-1;
cycles=[];
percent_error=[];
position=1; % reveals which loop is running

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
Force=zeros(N,N+1);
for k=1:length(y)
    for j=1:length(x)
        Force(k,j)=(sin(pi*(x(j)+pi)/(2*pi))*cos((pi/2)*(2*(y(k)+pi)/(2*pi)+1)));
    end
end

%% Gauss Seidel and Gauss Seidel with Relaxation Methods
finished=0; % tells the while loop when to stop

% The while loop below provokes the next loops to run with zero forcing function
while finished<=1
    lambda=1;
    
    % The while loop below switches the value of lambda from 1 to 1+"addon"
    % to provoke a second solution using the relaxation method.
    while lambda<=1+addon
        
        
        % The while loop below continually improves the solution by using
        % the newly solved values of the domain to solve the domain again,
        % until the error between iterations is less than 1 percent
        
        u=u_domain; % Sets or resets the domain boundary conditions and zeros out all others
        errorval=100; %  arbitrary initial error value
        iterations=0; % number of iterations for each solution, later stored as an element in the array "cycles"
        while errorval>1
            iterations=iterations+1;
            
            u_old=u; % This value is used to calculate error from iteration to iteration, and is also used in the relaxation equation
            
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
            errorval=mean(mean(error));
            u=lambda*u+(1-lambda)*u_old;
        end
        if position==1  % Stores solutions for later use
            u_gauss_seidel=u;
        elseif position==2
            u_gauss_seidel_relaxed=u;
        elseif position==3
            u_gauss_seidel_no_force=u;
        elseif position==4
            u_gauss_seidel_no_force_relaxed=u;
        else
            disp(error)
        end
        cycles(end+1)=iterations;
        percent_error(end+1)=errorval;
        lambda=lambda+addon;
        position=position+1;
    end
    lambda=1;
    Force=Force*0;
    finished=finished+1;
end


clear error u_old u_domain counter counter2 j k x y iterations addon Force N errorVal finished
end