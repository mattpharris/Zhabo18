Nonruled = csvread('OutputDutyCurves_kPP1_limit_256_ruled.TXT');
Ruled = csvread('OutputDutyCurves_kPP1_limit_256.TXT');

Labels = {};
i = 0;
for kPP1_ = [.25, .5, 1, 2, 4, 8, 16, 32, 64, 128, 256]
    for kCaMKII_ = [.05]
        i = 1;
        Labels(end+1) = {strcat('kCaMKII_{',num2str(kCaMKII_),'}...kPP1_{',num2str(kPP1_),'}...Ratio(',num2str(max(1,kCaMKII_/kPP1_)),':',num2str(max(1,kPP1_/kCaMKII_)),')')};
    end
end
LabelsRuled = strcat('Ruled... ',Labels);
LabelsNonruled = strcat('Nonruled... ',Labels);
Leg = [LabelsRuled LabelsNonruled];

%Plot The Duty Cycle for Ca2+ sweep
%plottingGroup = [1:7:length(Labels)];
plottingGroup = [1:length(Labels)];

c = num2cell(gray(length(plottingGroup)+1),2);
c = c(1:end-1);
figure;
hold on
subplot(2,1,1)
ruledPlot = plot(Ruled(:,1),Ruled(:,plottingGroup + 1),'LineWidth',2);
legend([LabelsRuled(plottingGroup)])
set(ruledPlot, {'color'}, c)
colormap(gray)
colorbar('Direction','reverse','Ticks',[0,1],'TickLabels',{'Low kCaMKII:kPP1 (1:160)','High kCaMII:kPP1 (16:1)'});
xlabel('[Ca2+] (uM)')
ylabel('Total Phosphorylated Subunits (uM)')

subplot(2,1,2)
nonruledPlot = plot(Nonruled(:,1),Nonruled(:,plottingGroup + 1),'LineWidth',2);
legend([LabelsNonruled(plottingGroup)])
set(nonruledPlot, {'color'}, c)
colormap(gray)
colorbar('Direction','reverse','Ticks',[0,1],'TickLabels',{'Low kCaMKII:kPP1 (1:160)','High kCaMII:kPP1 (16:1)'});
xlabel('[Ca2+] (uM)')
ylabel('Total Phosphorylated Subunits (uM)')