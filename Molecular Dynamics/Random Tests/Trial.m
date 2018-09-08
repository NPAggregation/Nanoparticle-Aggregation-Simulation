clear,clc
figure(1)
a=[0 0 1; 0 1 0; 1 0 0];
plot3(a(:,1),a(:,2),a(:,3),'or')
ylabel('y'), xlabel('x'),zlabel('z');
ylim([0,10]),xlim([0,10]),zlim([0,10]);
set(gca,'Color',[1 1 1])
for i=1:1:3
fprintf('Point coordinates %d %d %d\n',a(i,:))
end