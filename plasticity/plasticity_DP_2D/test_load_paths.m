level=3;

linewdth=2; 

elem_type='P1'; 
plasticity_DP_fem; 
zeta_hist_P1=zeta_hist(1:step); pressure_hist_P1=pressure_hist(1:step);

elem_type='P2'; 
plasticity_DP_fem; 
zeta_hist_P2=zeta_hist(1:step); pressure_hist_P2=pressure_hist(1:step);

elem_type='Q1'; 
plasticity_DP_fem; 
zeta_hist_Q1=zeta_hist(1:step); pressure_hist_Q1=pressure_hist(1:step);

elem_type='Q2'; 
plasticity_DP_fem; 
zeta_hist_Q2=zeta_hist(1:step); pressure_hist_Q2=pressure_hist(1:step);

figure
semilogx(zeta_hist_P1,pressure_hist_P1/c0,'-',...
         zeta_hist_P2,pressure_hist_P2/c0,'-',...
         zeta_hist_Q1,pressure_hist_Q1/c0,'-',...
         zeta_hist_Q2,pressure_hist_Q2/c0,'-','LineWidth',linewdth);
xlabel('displacement'); ylabel('normalized pressure');
legend('P1 elements','P2 elements','Q1 elements','Q2 elements','location','southeast') 
axis tight; enlarge_axis(0.05,0.05); 

if printout
      fig = gcf; fig.PaperPositionMode = 'auto'; fig_pos = fig.PaperPosition; fig.PaperSize = [fig_pos(3) fig_pos(4)];
      print('-painters','-dpdf',strcat('figures/fig_DP_2D_load_path_level_',num2str(level)))
end

