function show_grid_and_ground_truth(gridLF,source_idx);
gridLF.pos(gridLF.inside,:);
figure;ft_plot_mesh(ans)
hold all;ft_plot_mesh(ans(source_idx,:),'vertexcolor','r','vertexsize',100)
axis on
xlabel('x');
ylabel('y');
zlabel('z');
set(gca,'FontSize',16)
end