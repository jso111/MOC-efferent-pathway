
clear;
T = readtable('intensity ratio_summary1.xlsx');
T.Genotype=categorical(T.Genotype);
T.Genotype=reordercats(T.Genotype,{'WT','VGLUT3KO'});
T.Condition=categorical(T.Condition);
T.Condition=reordercats(T.Condition,{'OHClow','DClow','OHCmid','DCmid','OHChigh','DChigh'});

boxchart(T.Condition,T.Value,'GroupByColor',T.Genotype,'LineWidth',2);
set(gca,'FontSize',12)
legend('WT','VGLUT3^{-/-}','Location','northeast');
ylabel('Activity (re:RWM)','FontSize',14);
