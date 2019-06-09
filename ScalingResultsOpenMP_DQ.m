%% Dean Pakravan: 757389
% Qijie Li: 927249
% Scaling results in Open MP for the Shallow water equations
% Code - This was adapted from the ScalingResultsOpenMP.m file provided by
% Stephen Moore

% Barcoo - Delta_x=Delta_y=0.1, Delta_t = 0.2 (t_max = 100)
% Simulation took N/A seconds with 64 threads
% Simulation took N/A seconds with 32 threads 
% Simulation took 2277.585 seconds with 16 threads
% Simulation took 2392.15 seconds with 8 threads
% Simulation took 2316.59 seconds with 4 threads
% Simulation took 4935.65 seconds with 2 threads
% Simulation took 8486.25 seconds with 1 threads

%% Code
close all; clear all; clc;

saveFigures     =  1;
fontSize        = 30;
lineWidth       =  3;
markerSize      = 30;
paperSize       = [10, 10];
paperPosition   = [0, 0, 10, 10];

numThreads      = [1 2 4 8 16];
runTime         = [8486.25 4935.65 2316.59 2392.15 2277.585];

speedup         = runTime(1)./runTime;
ideal           = 2.^(0:length(numThreads)-1);
efficiency      = speedup./numThreads*numThreads(1);
log2NumThreads    = log2(numThreads);
log2Speedup     = log2(speedup);
log2Ideal       = log2(ideal);
efficiencyRange	= 0.4:0.1:1.1;
xtl             = {};
ytl             = {};

figure;
plot(log2NumThreads, log2Speedup, '.-', log2NumThreads, log2Ideal, '-', 'LineWidth', lineWidth, 'MarkerSize', markerSize);
grid on;
set(get(gca,'Xlabel'),  'String', 'Number of Threads', 'Fontsize', fontSize);
set(get(gca,'Ylabel'),  'String', 'Speedup',           'Fontsize', fontSize);
legend({'Observed', 'Ideal'}, 'Location', 'NorthWest');
axis([log2NumThreads(1), log2NumThreads(end), log2Ideal(1), log2Ideal(end)]);
set(gca, 'XTick', log2NumThreads);
set(gca, 'YTick', log2Ideal);
for j = 1:length(numThreads)
  xtl{j} = num2str(numThreads(j));
  ytl{j} = num2str(ideal(j));
end
set(gca, 'XTickLabel', xtl);
set(gca, 'YTickLabel', ytl);
set(gca, 'Fontsize',   fontSize);
set(gcf,'Color', [1,1,1],'PaperSize', paperSize, 'PaperPosition', paperPosition);
if saveFigures
	saveas(gcf,['ShallowWaterEquationsOpenMPSpeedup.pdf']);
end

figure;
plot(log2NumThreads, efficiency, '.-', log2NumThreads, ones(length(numThreads),1), '-', 'LineWidth', lineWidth, 'MarkerSize', markerSize);
grid on;
set(get(gca,'Xlabel'),  'String', 'Number of Nodes',     'Fontsize', fontSize);
set(get(gca,'Ylabel'),  'String', 'Parallel Efficiency', 'Fontsize', fontSize);
legend({'Observed', 'Ideal'}, 'Location', 'NorthWest');
axis([log2NumThreads(1), log2NumThreads(end), efficiencyRange(1), efficiencyRange(end)]);
set(gca, 'XTick', log2NumThreads);
set(gca, 'YTick', efficiencyRange);
for j = 1:length(numThreads)
  xtl{j} = num2str(numThreads(j));
end
for j = 1:length(efficiencyRange)
    ytl{j} = num2str(efficiencyRange(j), '%0.1f');
end
set(gca, 'XTickLabel', xtl);
set(gca, 'YTickLabel', ytl);
set(gca, 'Fontsize',   fontSize);
set(gcf,'Color', [1,1,1],'PaperSize', paperSize, 'PaperPosition', paperPosition);
if saveFigures
	saveas(gcf,['ShallowWaterEquationsOpenMPParallelEfficiency.pdf']);
end