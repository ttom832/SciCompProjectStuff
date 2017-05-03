clear
clc
Lambda=1.4; % verified with experimentation to be best value
y=[10]; % geometric multiplier which increases number of nodes in domain
Nodes=y*10; % Nodes along edge of sqaure domain
analysisNode1=2;
analysisNode2=9;

for p=1:length(Nodes)
[h,...
    X,...
    Y,...
    cycles,...
    percent_error,...
    u_gauss_seidel,...
    u_gauss_seidel_relaxed,...
    u_gauss_seidel_no_force,...
    u_gauss_seidel_no_force_relaxed]=project(Nodes(p),Lambda);


Cycles(p,:)=cycles;
Error(p,:)=percent_error;
H(p,1)=h;

% (-1.745 , 1.745) position
upper_left_position(p)=u_gauss_seidel(analysisNode1*y(p),analysisNode1*y(p)+1);
% (-1.745 , 1.745) position
upper_right_position(p)=u_gauss_seidel(analysisNode1*y(p),analysisNode2*y(p)+1);
% (-1.745 , 1.745) position
lower_left_position(p)=u_gauss_seidel(analysisNode2*y(p),analysisNode1*y(p)+1);
% (-1.745 , 1.745) position
lower_right_position(p)=u_gauss_seidel(analysisNode2*y(p),analysisNode2*y(p)+1);
end

upper_left_position
upper_right_position
lower_left_position
lower_right_position
disp('Step Size')
disp(H)
disp('Gauss Seidel , Gauss Seidel with Relaxation , Gauss Seidel no Forcing Function , Gauss Seidel with Relaxation no Forcing Function')
Cycles
Error
for k=1:length(Nodes)
percent_improvement(k,1:2)=[(Cycles(k,1)-Cycles(k,2))/(Cycles(k,1))*100,(Cycles(k,3)-Cycles(k,4))/(Cycles(k,3))*100];
end
percent_improvement


plot(Nodes,upper_right_position,Nodes,upper_left_position,Nodes,lower_left_position,Nodes,lower_right_position),legend('Position A','Position B','Position C','Position D'),grid,title('Grid Convergence'),xlabel('Nodes along Edge'),ylabel('Value at Position')
plot(Nodes,Cycles(:,1),Nodes,Cycles(:,2)),grid,legend('Gauss Seidel Cycles','Gauss Seidel Relaxation Cycles','Location','SouthEast'),xlabel('Nodes Along Edge'),ylabel('Cycles'),title('Cylces vs. Domain Density')
figure(1)
mesh(X,Y,u_gauss_seidel),xlabel('x'),ylabel('y'),zlabel('u'),title('Gauss Seidel')

figure(2)
mesh(X,Y,u_gauss_seidel_relaxed),xlabel('x'),ylabel('y'),zlabel('u'),title('Gauss Seidel With Relaxation')

figure(3)
mesh(X,Y,u_gauss_seidel_no_force),xlabel('x'),ylabel('y'),zlabel('u'),title('Gauss Seidel No Forcing Function')

figure(4)
mesh(X,Y,u_gauss_seidel_no_force_relaxed),xlabel('x'),ylabel('y'),zlabel('u'),title('Gauss Seidel No Forcing Function With Relaxation')

figure(5)
mesh(X,Y,u_gauss_seidel_relaxed-u_gauss_seidel),xlabel('x'),ylabel('y'),zlabel('u'),title('Difference in Gauss Seidel and Gauss Seidel Relaxed')














