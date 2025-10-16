function [] = generateAnimation(P,dt, strideLength)
figure()

Px = P(:,[1 3 5 7 9]);
Py = P(:,[2 4 6 8 10]); 

H=line(Px(1,:), Py(1,:)); hold on
axis equal
H=handle(H);
xlim([-0.6 3.6])
ylim([-0.1 0.8])
xlabel('(m)')
ylabel('(m)')

for i = 1:7
for j=1:size(P,1)
    H.XData=Px(j,:);
    H.YData=Py(j,:);
    pause(dt)
end
Px = Px + strideLength;

end
