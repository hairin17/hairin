% Shen network matrix plot
% input: loooverlap total matrix of power 268 aligned
% output:  # of edges by each network, Matrix figure

function [shenmat] = shensoc_matplot(overlap)

% shenlab = xlsread('C:\Users\Seyul\Dropbox\Researches\Researches\A_journal_Neuromarker\cpm\shen_268_parcellation_networklabels_SORTED.csv');
% shenlab = xlsread('C:\Program Files\MATLAB\R2015b\toolbox\spm12\spm12\toolbox\conn\rois\shen_268_parcellation_networklabels_SORTED.csv');

% shennet=shenlab(:,5);
% shennet = [2,4,3,2,3,3,2,2,2,1,4,1,3,2,4,1,2,4,2,4,2,2,5,5,5,5,5,4,4,2,2,4,5,5,5,4,5,5,5,5,8,6,8,4,5,5,2,2,3,3,5,1,1,1,2,1,1,5,8,5,5,5,5,1,1,8,8,6,8,2,8,6,8,8,6,7,6,7,6,6,7,6,4,5,3,3,6,4,5,3,4,5,4,4,4,3,5,6,4,7,4,7,4,4,4,4,4,4,5,4,2,2,4,4,3,2,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,4,4,1,3,2,1,3,2,2,4,1,4,2,1,1,1,1,4,1,2,4,1,2,5,5,5,5,1,5,2,1,5,5,5,4,5,5,5,5,5,8,6,8,4,5,5,5,2,1,2,1,1,1,5,5,1,5,1,2,1,5,2,5,6,2,8,8,5,3,8,6,8,6,6,8,8,6,7,7,7,6,6,4,5,1,4,4,3,3,4,3,4,3,5,4,4,4,4,4,4,5,4,4,4,3,8,7,2,4,4,4,2,2,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4]';
% socog=[1 2 3 5 10 16 19	49	50	53	55	64	65	90	99	137	138	140	145	148	183	185	187	190	222	225	231];%	242];  % 242 cereb
% shennet(socog) = 9;
shennet = [9,9,9,2,9,3,2,2,2,9,4,1,3,2,4,9,2,4,9,4,2,2,5,5,5,5,5,4,4,2,2,4,5,5,5,4,5,5,5,5,8,6,8,4,5,5,2,9,9,3,5,1,1,1,9,1,1,5,8,5,5,5,5,9,9,8,8,6,8,2,8,6,8,8,6,7,6,7,6,6,7,6,4,5,3,3,6,4,5,9,4,5,4,4,4,3,5,6,4,7,4,7,4,4,4,4,4,4,5,4,2,2,4,4,3,2,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,3,4,4,9,9,2,9,3,2,2,4,9,4,2,9,1,1,1,4,1,2,4,1,2,5,5,5,5,1,5,2,1,5,5,5,4,5,5,5,5,5,8,6,8,4,5,5,5,2,9,2,9,1,9,5,5,9,5,1,2,1,5,2,5,6,2,8,8,5,3,8,6,8,6,6,8,8,6,7,7,7,6,6,4,5,1,4,4,9,3,4,9,4,3,5,4,4,9,4,4,4,5,4,4,4,3,8,7,9,4,4,4,2,2,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4];

a{1} = find(shennet == 1);% MedialFrontal
a{2} = find(shennet == 2); % FrontoParietal
a{3} = find(shennet == 3);% Defaultmode
a{4} = find(shennet == 4);% SubcortCereb
a{5} = find(shennet == 5);% Motor
a{6} = find(shennet == 6);% Visual1
a{7} = find(shennet == 7);% Visual2
a{8} = find(shennet == 8);% VisualAssoc
a{9} = find(shennet == 9); % social cognition


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

% shenmat = tril(cell2);   %%%% Select sum
shenmat = tril(cell3);   %%%% Select mean

% figure(1);imagesc(cell2);colormap(autumn)  % pos
% figure(2);imagesc(cell2);colormap(winter)  % neg

imagesc(shenmat);            % Create a colored plot of the matrix values
colormap(flipud(gray));  % Change the colormap to gray (so higher values are
                         %   black and lower values are white)

textStrings = num2str(shenmat(:), '%0.2f');       % Create strings from the matrix values
% textStrings = num2str(shenmat(:));       % Create strings from the matrix values

textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
[x, y] = meshgrid(1:no_nets);  % Create x and y coordinates for the strings
hStrings = text(x(:), y(:), textStrings(:), ...  % Plot the strings
                'HorizontalAlignment', 'center');
% hStrings = tril(hStrings);
midValue = mean(get(gca, 'CLim'));  % Get the middle value of the color range
textColors = repmat(shenmat(:) > midValue, 1, 3);  % Choose white or black for the
                                               %   text color of the strings so
                                               %   they can be easily seen over
                                               %   the background color
set(hStrings, {'Color'}, num2cell(textColors, 2));  % Change the text colors
set(gca, 'XTick', 1:no_nets, ...                             % Change the axes tick marks
         'XTickLabel', {'MedF', 'FP', 'DMN', 'Sub', 'Mo', 'MVis', 'LVis', 'VisA', 'Social'}, ...  %   and tick labels
         'YTick', 1:no_nets, ...
         'YTickLabel', {'MedF', 'FP', 'DMN', 'Sub', 'Mo', 'MVis', 'LVis', 'VisA', 'Social'}, ...
         'TickLength', [0 0]);
     
     
end


