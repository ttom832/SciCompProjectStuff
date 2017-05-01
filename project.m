clear
clc

%% Gauss Seidel
N=100;
h=2*pi/(N-1);
x=linspace(-pi-h,pi,N+1);
y=linspace(pi,-pi,N);
u=zeros(N,N+1);
[X,Y]=meshgrid(x,y);
% Top B.C.'s
% f_ax
for k=2:N+1
    u(1,k)=x(k)*(x(k)+pi)^2;
end

% Bottom B.C.'s
% g_ax
for k=2:N+1
    u(N,k)=(x(k)+pi)^2*cos(-x(k));
end

% Right B.C.'s
% H_a
for k=1:N
    u(k,N+1)=((2*pi)^2*cos(-pi))+...
        ((y(k)+pi)/(2*pi))*...
        (pi*(2*pi)^2-...
        ((2*pi)^2*cos(-pi)));
end


errorval=100;
iterations=0;
while errorval>1
    iterations=iterations+1;
    u_old=u;
    % Circling
    counter=2;
    while counter<=N/2
        
        % Top Pyramid
        for j=      N-counter+2   :  -1 :   counter+1
            u( counter , j ) = ((sin(pi*(x(j)+pi)/(2*pi))*cos((pi/2)*(2*(y(counter)+pi)/(2*pi)+1)))*h^2+...
                u(counter-1,j)+u(counter+1,j)+u(counter,j-1)+u(counter,j+1))*0.25;
        end
        
        % Bottom Pyramid
        for j=      N-counter+2  :  -1 :    counter+1
            u( N-counter+1 , j ) = ((sin(pi*(x(j)+pi)/(2*pi))*cos((pi/2)*(2*(y(N-counter+1)+pi)/(2*pi)+1)))*h^2+...
                u(N-counter+2,j)+u(N-counter,j)+u(N-counter+1,j-1)+u(N-counter+1,j+1))*0.25;
        end
        
        % Right Pyramid
        for k=       counter+1      :        N-counter
            u( k , N-counter+2 ) = ((sin(pi*(x(N-counter+2)+pi)/(2*pi))*cos((pi/2)*(2*(y(k)+pi)/(2*pi)+1)))*h^2+...
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
                u( k , j) = ((sin(pi*(x(j)+pi)/(2*pi))*cos((pi/2)*(2*(y(k)+pi)/(2*pi)+1)))*h^2+...
                    u(k-1,j)+u(k+1,j)+u(k,j-1)+u(k,j+1))*0.25;
            end
            
        end
        
    else % odd number of nodes
        for j=       floor(N/2)+2 :   -1 :  3
            counter2=counter2+1;
            for k =  floor(N/2)-counter2+2   :    floor(N/2)+counter2
                u( k , j) = ((sin(pi*(x(j)+pi)/(2*pi))*cos((pi/2)*(2*(y(k)+pi)/(2*pi)+1)))*h^2+...
                    u(k-1,j)+u(k+1,j)+u(k,j-1)+u(k,j+1))*0.25;
            end
            
        end
    end
    
    % ghost nodes
    u(:,1)=u(:,3);
    
    % left boundary nodes
    for k=2:N-1
        u(k,2)=((sin(pi*(x(2)+pi)/(2*pi))*cos((pi/2)*(2*(y(k)+pi)/(2*pi)+1)))*h^2+...
            u(k-1,2)+u(k+1,2)+u(k,1)+u(k,3))*0.25;
    end
    
    % remaining 2 nodes at nueman boundary
    u(1,2)=((sin(pi*(x(2)+pi)/(2*pi))*cos((pi/2)*(2*(y(1)+pi)/(2*pi)+1)))*h^2+...
        u(2,2)+u(1,1)+u(1,3))*0.25;
    u(N,2)=((sin(pi*(x(2)+pi)/(2*pi))*cos((pi/2)*(2*(y(N)+pi)/(2*pi)+1)))*h^2+...
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
end
iterations
errorval
figure(1)
mesh(X,Y,u),xlabel('x'),ylabel('y'),zlabel('u'),title('Gauss Seidel Method')

u_rel=zeros(N,N+1);


%% Gauss-Seidel with Relaxation

% Top B.C.'s
for k=2:N+1
    u_rel(1,k)=x(k)*(x(k)+pi)^2;
end

% Bottom B.C.'s
for k=2:N+1
    u_rel(N,k)=(x(k)+pi)^2*cos(pi*x(k)/-pi);
end

% Right B.C.'s
for k=1:N
    u_rel(k,N+1)=((pi+pi)^2*cos(pi*pi/-pi))+...
        ((y(k)+pi)/(pi+pi))*...
        (pi*(pi+pi)^2-...
        ((pi+pi)^2*cos(pi*pi/-pi)));
end

errorval2=100;
iterations2=0;
lambda=1.4;
while errorval2>1
    iterations2=iterations2+1;
    u_old2=u_rel;
    % Circling
    counter=2;
    while counter<=N/2
        
        % Top Pyramid
        for j=      N-counter+2   :  -1 :   counter+1
            u_rel( counter , j ) = ((sin(pi*(x(j)+pi)/(2*pi))*cos((pi/2)*(2*(y(counter)+pi)/(2*pi)+1)))*h^2+...
                u_rel(counter-1,j)+u_rel(counter+1,j)+u_rel(counter,j-1)+u_rel(counter,j+1))*0.25;
        end
        
        
        % Bottom Pyramid
        for j=      N-counter+2  :  -1 :    counter+1
            u_rel( N-counter+1 , j ) = ((sin(pi*(x(j)+pi)/(2*pi))*cos((pi/2)*(2*(y(N-counter+1)+pi)/(2*pi)+1)))*h^2+...
                u_rel(N-counter+2,j)+u_rel(N-counter,j)+u_rel(N-counter+1,j-1)+u_rel(N-counter+1,j+1))*0.25;
        end
        
        % Right Pyramid
        for k=       counter+1      :        N-counter
            u_rel( k , N-counter+2 ) = ((sin(pi*(x(N-counter+2)+pi)/(2*pi))*cos((pi/2)*(2*(y(k)+pi)/(2*pi)+1)))*h^2+...
                u_rel(k-1,N-counter+2)+u_rel(k+1,N-counter+2)+u_rel(k,N-counter+1)+u_rel(k,N-counter+3))*0.25;
        end
        
        
        counter=counter+1;
        
    end
    
    % Left Pyramid
    counter2=0;
    if mod(N,2)==0 % even number of nodes
        for j=       floor(N/2)+1 :   -1 :  3
            counter2=counter2+1;
            for k =  floor(N/2)-counter2+2    :    floor(N/2)+counter2-1
                u_rel( k , j) = ((sin(pi*(x(j)+pi)/(2*pi))*cos((pi/2)*(2*(y(k)+pi)/(2*pi)+1)))*h^2+...
                    u_rel(k-1,j)+u_rel(k+1,j)+u_rel(k,j-1)+u_rel(k,j+1))*0.25;
            end
            
        end
        
    else % odd number of nodes
        for j=       floor(N/2)+2 :   -1 :  3
            counter2=counter2+1;
            for k =  floor(N/2)-counter2+2   :    floor(N/2)+counter2
                u_rel( k , j) = ((sin(pi*(x(j)+pi)/(2*pi))*cos((pi/2)*(2*(y(k)+pi)/(2*pi)+1)))*h^2+...
                    u_rel(k-1,j)+u_rel(k+1,j)+u_rel(k,j-1)+u_rel(k,j+1))*0.25;
            end
            
        end
    end
    
    % ghost nodes
    u_rel(:,1)=u_rel(:,3);
    
    % left boundary nodes
    for k=2:N-1
        u_rel(k,2)=((sin(pi*(x(2)+pi)/(2*pi))*cos((pi/2)*(2*(y(k)+pi)/(2*pi)+1)))*h^2+...
            u_rel(k-1,2)+u_rel(k+1,2)+u_rel(k,1)+u_rel(k,3))*0.25;
    end
    
    % remaining 2 nodes at nueman boundary
    u_rel(1,2)=((sin(pi*(x(2)+pi)/(2*pi))*cos((pi/2)*(2*(y(1)+pi)/(2*pi)+1)))*h^2+...
        u_rel(2,2)+u_rel(1,1)+u_rel(1,3))*0.25;
    u_rel(N,2)=((sin(pi*(x(2)+pi)/(2*pi))*cos((pi/2)*(2*(y(N)+pi)/(2*pi)+1)))*h^2+...
        u_rel(N-1,2)+u_rel(N,1)+u_rel(N,3))*0.25;
    
    
    for k=1:N
        for j=1:N+1
            error2(k,j)=abs((u_rel(k,j)-u_old2(k,j))/u_rel(k,j))*100;
        end
    end
    % average error over domain
    errorval2=mean(mean(error2));
    u_rel=lambda*u_rel+(1-lambda)*u_old2;
end
iterations2
errorval2
figure(2)
mesh(X,Y,u_rel),xlabel('x'),ylabel('y'),zlabel('u'),title('Gauss Seidel with Relaxation')

%% Gauss Siedel No Forcing Function
uNoF=zeros(N,N+1);
% Top B.C.'s
for k=2:N+1
    uNoF(1,k)=x(k)*(x(k)+pi)^2;
end

% Bottom B.C.'s
for k=2:N+1
    uNoF(N,k)=(x(k)+pi)^2*cos(pi*x(k)/-pi);
end

% Right B.C.'s
for k=1:N
    uNoF(k,N+1)=((pi+pi)^2*cos(pi*pi/-pi))+...
        ((y(k)+pi)/(pi+pi))*...
        (pi*(pi+pi)^2-...
        ((pi+pi)^2*cos(pi*pi/-pi)));
end

errorval3=100;
iterations3=0;
while errorval3>1
    iterations3=iterations3+1;
    u_oldNoF=uNoF;
    % Circling
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
    errorval3=mean(mean(error3));
end
iterations3
errorval3
figure(3)
mesh(X,Y,uNoF),xlabel('x'),ylabel('y'),zlabel('u'),title('Gauss Seidel with no Forcing Function')
