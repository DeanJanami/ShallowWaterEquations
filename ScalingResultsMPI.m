% Avoca - myN_x=myN_y=1000, N_t=100,000

% Simulation took 70.1671 seconds with 4 processes
% Simulation took 70.3472 seconds with 9 processes
% Simulation took 75.3253 seconds with 16 processes
% Simulation took 75.3509 seconds with 25 processes
% Simulation took 75.3962 seconds with 36 processes
% Simulation took 75.3944 seconds with 49 processes
% Simulation took 75.3441 seconds with 64 processes
% Simulation took 75.3686 seconds with 81 processes
% Simulation took 75.3535 seconds with 100 processes

close all; clear all; clc;

saveFigures     =  1;
fontSize        = 30;
lineWidth       =  3;
markerSize      = 30;
paperSize       = [10, 10];
paperPosition   = [0, 0, 10, 10];

numProcs    	= [4 9 16 25 36 49 64 81 100];
runTime         = [70.1671 70.3472 75.3253 75.3509 75.3962 75.3944 75.3441 75.3686 75.3535];

figure;
plot(numProcs, runTime, '.-', 'LineWidth', lineWidth, 'MarkerSize', markerSize);
grid on;
set(get(gca,'Xlabel'),  'String', 'Number of Processes', 'Fontsize', fontSize);
set(get(gca,'Ylabel'),  'String', 'Run Time [s]',        'Fontsize', fontSize);
legend({'Observed', 'Ideal'}, 'Location', 'NorthWest');
% axis([log2NumThreads(1), log2NumThreads(end), log2Ideal(1), log2Ideal(end)]);
% set(gca, 'XTick', log2NumThreads);
% set(gca, 'YTick', log2Ideal);
% for j = 1:length(numThreads)
%   xtl{j} = num2str(numThreads(j));
%   ytl{j} = num2str(ideal(j));
% end
% set(gca, 'XTickLabel', xtl);
% set(gca, 'YTickLabel', ytl);
set(gca, 'Fontsize',   fontSize);
set(gcf,'Color', [1,1,1],'PaperSize', paperSize, 'PaperPosition', paperPosition);
if saveFigures
	saveas(gcf,['ShallowWaterEquationsMPI Run Time.pdf']);
end