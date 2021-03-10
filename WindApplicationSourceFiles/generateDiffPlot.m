%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DESCRIPTION:
%   Script to generate power curve difference plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for pair = 1:6
    p = [];
    hold off;
    p(1) = plot(xtest,muDiff(:,pair),'LineWidth',2,'Color','blue');
    hold on;
    p(2) = plot(xtest,band(:,pair),'--','LineWidth',2,'Color','red');
    p(3) = plot(xtest,-band(:,pair),'--','LineWidth',2,'Color','red');
    xlabel('Wind speed (m/s)','FontSize',18);
    xlim([2,16]);
    ylim([-20,20]);
    legend([p(1),p(2)],'Difference','90% Confidence band','Location','northeast'); 
    set(gca,'FontSize',18);
    set(gcf,'PaperSize', [6 5]);
    hold off;
    if pair == 1
        plotmain = "Year 2 vs Year 1";
        pairname = "Figure4a";
    elseif pair == 2
        plotmain = "Year 3 vs Year 1";
        pairname = "Figure4b";
    elseif pair == 3
        plotmain = "Year 4 vs Year 1";
        pairname = "Figure4c";
    elseif pair == 4
        plotmain = "Year 3 vs Year 2";
        pairname = "Figure4d";
    elseif pair == 5
        plotmain = "Year 4 vs Year 2";
        pairname = "Figure4e";
    elseif pair == 6
        plotmain = "Year 4 vs Year 3";
        pairname = "Figure4f";
    end
    
    title(plotmain);
    plotname = strcat(pairname,".pdf");
    print('-painters',plotname,'-dpdf');
end
    
    
    
    