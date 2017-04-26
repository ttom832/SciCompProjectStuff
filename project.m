clear
clc
%%
N=10;
h=2*pi/(N-1);
x=linspace(-pi-h,pi,N+1);
y=linspace(pi,-pi,N);
u=zeros(N,N+1);
[X,Y]=meshgrid(x,y);


% Top B.C.'s
for k=2:N+1
    u(1,k)=x(k)*(x(k)+pi)^2;
end

% Bottom B.C.'s
for k=2:N+1
    u(N,k)=(x(k)+pi)^2*cos(pi*x(k)/-pi);
end

% Right B.C.'s
for k=1:N
    u(k,N+1)=((pi+pi)^2*cos(pi*pi/-pi))+...
        ((y(k)+pi)/(pi+pi))*...
        (pi*(pi+pi)^2-...
        ((pi+pi)^2*cos(pi*pi/-pi)));
end


iterations=0;
while iterations<5
    iterations=iterations+1
    u_old=u;
    % Circling
    counter=2;
    while counter<=N/2

        % Top Pyramid
        for j=      N-counter+2   :  -1 :   counter+1          
        u( counter , j ) = ((sin(pi*(x(j)+pi)/(2*pi))*cos((pi/2)*(2*(y(counter)+pi)/(2*pi)+1)))*h^2+...
             u(counter-1,j)+u(counter+1,j)+u(counter,j-1)+u(counter,j+1))/4;
        end


        % Bottom Pyramid
        for j=      N-counter+2  :  -1 :    counter+1        
            u( N-counter+1 , j ) = ((sin(pi*(x(j)+pi)/(2*pi))*cos((pi/2)*(2*(y(N-counter+1)+pi)/(2*pi)+1)))*h^2+...
             u(N-counter+2,j)+u(N-counter,j)+u(N-counter+1,j-1)+u(N-counter+1,j+1))/4;
        end

        % Right Pyramid
        for k=       counter+1      :        N-counter
            u( k , N-counter+2 ) = ((sin(pi*(x(N-counter+2)+pi)/(2*pi))*cos((pi/2)*(2*(y(k)+pi)/(2*pi)+1)))*h^2+...
             u(k-1,N-counter+2)+u(k+1,N-counter+2)+u(k,N-counter+1)+u(k,N-counter+3))/4;
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
               u(k-1,j)+u(k+1,j)+u(k,j-1)+u(k,j+1))/4;
        end

    end

    else % odd number of nodes
    for j=       floor(N/2)+2 :   -1 :  3
        counter2=counter2+1;
        for k =  floor(N/2)-counter2+2   :    floor(N/2)+counter2
            u( k , j) = ((sin(pi*(x(j)+pi)/(2*pi))*cos((pi/2)*(2*(y(k)+pi)/(2*pi)+1)))*h^2+...
               u(k-1,j)+u(k+1,j)+u(k,j-1)+u(k,j+1))/4;
        end

    end
    end

    % ghost nodes
    u(:,1)=u(:,3)

    % left boundary nodes
    for k=2:N-1
        u(k,2)=((sin(pi*(x(2)+pi)/(2*pi))*cos((pi/2)*(2*(y(k)+pi)/(2*pi)+1)))*h^2+...
               u(k-1,2)+u(k+1,2)+u(k,1)+u(k,3))/4;
    end

    % remaining 2 nodes at nueman boundary
    u(1,2)=((sin(pi*(x(2)+pi)/(2*pi))*cos((pi/2)*(2*(y(1)+pi)/(2*pi)+1)))*h^2+...
           u(2,2)+u(1,1)+u(1,3))/4;
    u(N,2)=((sin(pi*(x(2)+pi)/(2*pi))*cos((pi/2)*(2*(y(N)+pi)/(2*pi)+1)))*h^2+...
           u(N-1,2)+u(N,1)+u(N,3))/4;

       
       for k=1:N
           for j=1:N+1
   error(k,j)=abs((u(k,j)-u_old(k,j))/u(k,j))*100;
           end
       end
       error
end

mesh(u)


























