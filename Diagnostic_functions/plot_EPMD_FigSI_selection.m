t_opt = linspace(-2,36,77);
w=6;

t=linspace(-2,36,1000)';

s = exp(-((t-t_opt)./w).^2);

plot(t,s,'LineWidth',1)
set(gca,'FontSize',18)
axis([-2 36 0 1])
xlabel('Environmental Temperature (\circC)','FontSize',18)
ylabel('Selection Coefficient','FontSize',18)


sname=[pathname '../Figures/Figure_SI_selection.png'];
set(gcf,'Color','w')
export_fig(sname,'-r300')