function AnimatePlot(N, Particle, side)
pos = zeros(N, 3);

for i = 1:3
    pos(:, i) = GetVectorProps(Particle, N, i);
end

figure
plot3(pos(:,1), pos(:,2), pos(:,3), 'o', 'MarkerFaceColor', 'k');
xlabel('x'), ylabel('y'), zlabel('z'), title('Coordinates of Particles');
xlim([0 ceil(side)]), ylim([0 ceil(side)]), zlim([0 ceil(side)]);
grid on
drawnow

end