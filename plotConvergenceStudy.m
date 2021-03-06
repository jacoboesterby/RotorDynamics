function plotConvergenceStudy(folder,Exp,setup)

files = dir(fullfile(folder, '*.mat'));
files = {files.name};
value = zeros(1,length(files));
for k=1:length(files)
     B = regexp(files{k},'\d*','Match');
     value(k) = str2num(B{1});
end
figure
hold on
color = ['r','g','m','b','c'];
for i=1:length(value)
    load(strcat(folder,'/',files{i}))
    for k=1:length(eigFreq)
    plot(value(i),eigFreq(k),'*','color',color(k),'markersize',18,'HandleVisibility','off')
        if i==length(value)
        text(value(i)-2,eigFreq(k)+150,sprintf('e= %.1f %%',((Exp(k)-eigFreq(k))/Exp(k)*100)),'Fontsize',26,'fontname','times','color',color(k))
        end
    end

end
leg={};
for i=1:length(Exp)
    plot([0,max(value)+4],[Exp(i),Exp(i)],'color',color(i),'linewidth',3);
    leg{i} = sprintf('$f_{%d}$ = %.0f Hz',i,Exp(i));
end
legend(leg,'interpreter','latex','location','best');
xlabel('Number of elements [-]','interpreter','latex')
ylabel('Eigen frequency f [Hz]','interpreter','latex')
title(strcat('Convergence study - ',setup))
set(gca,'fontsize',28);
set(gcf, 'Position',  [0, 0, 1280, 800]);
saveas(gcf,strcat(folder,'/','Convergence study.png'))
end
