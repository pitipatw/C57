A = 10*rand(2,2);

b = 5*rand(2,1); 

x = rand(2,1,10000);
cones = zeros(2,1,size(x,3));
hyperplane = [] 
close all 
figure(1)
hold on
for i = 1:size(x,3) 
    xi = x(:,:,i);
    cones(:,1,i) = A*xi;
end

q = quiver(cones(1,1,:), cones(2,1,:))
q.ShowArrowHead = 'off';
q.Marker = 'none';


quiver(b(1),b(2), 'r', 'LineWidth',10)

%find p that pTb < 0 
p1 = rand(1) ;
p2 = (-rand(1) - p1*b(1))/b(2);
p = [p1 ;p2]
pa = p'*A
quiver(p1/norm(p),p2/norm(p), 'g')
quiver(pa(1), pa(2), 'k')
axis equal
