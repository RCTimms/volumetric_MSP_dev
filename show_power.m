function show_power(X,gridLF,source_idx,log_it)
power=sum(X.^2,2);
title_string='';
if log_it==1
power=log(power);
title_string='Log of ';
end
figure;
scatter3(gridLF.pos(gridLF.inside,1),gridLF.pos(gridLF.inside,2),gridLF.pos(gridLF.inside,3),50,power,'filled')
hold all;

gridLF.pos(gridLF.inside,:);

hold all;ft_plot_mesh(ans(source_idx,:),'vertexcolor','r')
colorbar
colormap jet
title([title_string,'Reconstructed Power']);
legend('','','Ground Truth Location');
end