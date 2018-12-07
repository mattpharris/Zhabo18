Ruled = csvread('PP1-CaMKII_Ruled_12-06_06-00.TXT');
Nonruled = csvread('PP1-CaMKII_Nonruled_12-06_06-00.TXT');

Labels = {}
i = 0
for kPP1_ = [.25, .5, 1, 2, 4, 8, 16]
    for kCaMKII_ = [.05, .1, .25, .5, 1, 2, 4]
        i = 1;
        Labels(end+1) = {strcat('kPP1_{',num2str(kPP1_),'}...kCaMKII_{',num2str(kCaMKII_),'}')};
    end
end
LabelsRuled = strcat('Ruled... ',Labels);
LabelsNonruled = strcat('Nonruled... ',Labels);
Leg = [LabelsRuled LabelsNonruled];

%Plot The Duty Cycle for Ca2+ sweep
plottingGroup = [4:14:length(Labels)];

figure
hold on
plot(Ruled(:,1),Ruled(:,plottingGroup + 1),'r');
plot(Nonruled(:,1),Nonruled(:,plottingGroup + 1),'b');
legend([LabelsRuled(plottingGroup) LabelsNonruled(plottingGroup)])
xlabel('[Ca2+] (uM)')
ylabel('Total Phosphorylated Subunits (uM)')