% Power network matrix plot
% input: loooverlap total matrix of power 256 aligned
% output: [sum, mean] # of edges by each network

function [powermat] = powedge_matplot(overlap)

% overlap = readtable('C:\Users\Seyul\Downloads\Telegram Desktop\overlapping mats\01_task-task\overlap_pos_100per.txt');


% a: creat network mask cells
a{1} = [1:35]; % 1:30 Sensory 31:35 mouth
a{2} = [36:49]; % 36: 49 con
a{3} = [50:62];% 51: 62 audi
a{4} = [63:125];% 63:120 dmn 121:125 memory
a{5} = [126:156];% 126:156 vis
a{6} = [157:181];% 157:181 FP
a{7} = [182:199];% 182:199 salience
a{8} = [200:212];% 200:212 subcort
a{9} = [213: 232];% 213: 232  ventral / dorsal attention
%exclude cerebell / uncertain

no_nets = length(a);

% assign each network mask into m by n cell
cell = {};
for m = 1:no_nets
    for n =1:no_nets
        cell{m,n} = overlap(a{m}, a{n});
    end
end

% sum within each cell
cell2=zeros(no_nets,no_nets);
for m = 1:no_nets
    for n = 1:no_nets
        cell2(m,n)=sum(sum(cell{m,n}));
    end
end

% mean within each cell
cell3=zeros(no_nets,no_nets);
for m = 1:no_nets
    for n = 1:no_nets
        cell3(m,n)=mean(mean(cell{m,n}))*100;
%         cell3(m,n)=sum(sum(cell{m,n}))/(size(cell{m,n},1)*size(cell{m,n},2));
    end
end

powermat = tril(cell2);   %%%% Select sum
% powermat = tril(cell3);   %%%% Select mean

imagesc(powermat);            % Create a colored plot of the matrix values
colormap(flipud(bone));  % Change the colormap to gray (so higher values are
                         %   black and lower values are white)

textStrings = num2str(powermat(:), '%0.2f');       % Create strings from the matrix values  '%0.2f'
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
[x, y] = meshgrid(1:no_nets);  % Create x and y coordinates for the strings
hStrings = text(x(:), y(:), textStrings(:), ...  % Plot the strings
                'HorizontalAlignment', 'center');
midValue = mean(get(gca, 'CLim'));  % Get the middle value of the color range
textColors = repmat(powermat(:) > midValue, 1, 3);  % Choose white or black for the
                                               %   text color of the strings so
                                               %   they can be easily seen over
                                               %   the background color
set(hStrings, {'Color'}, num2cell(textColors, 2));  % Change the text colors
set(gca, 'XTick', 1:no_nets, ...                             % Change the axes tick marks
         'XTickLabel', {'Sen', 'CO', 'Aud', 'DMN', 'Vis', 'FP', 'Sal', 'Sub', 'Att'}, ...  %   and tick labels
         'YTick', 1:no_nets, ...
         'YTickLabel', {'Sen', 'CO', 'Aud', 'DMN', 'Vis', 'FP', 'Sal', 'Sub', 'Att'}, ...
         'TickLength', [0 0]);
     
     
end


