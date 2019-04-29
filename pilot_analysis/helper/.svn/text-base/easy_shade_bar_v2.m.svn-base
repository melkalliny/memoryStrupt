function easy_shade_bar_v2(varargin)
% TS 2016
% shadedErrorbar adapted from JW
% easy_shade_bar(data,more_data,...)
% data rows should be examples that will be averaged across
% eg. epochs x time

alpha = .4;
n_in = nargin;
if ndims(varargin{1}) == 3
    foo =varargin{1};
    n_grps = size(foo,1);
    for i = 1:n_grps
        varargin{i} = squeeze(foo(i,:,:));
    end
    n_in = n_grps;
end
n_plots = n_in;

inputColor = varargin{2};

%color_wheel = {'r','b',[.6 0 1], 'c','m','g','k','y'};
c1 = [.7 0 0];
c2 = [.9 .7 .1];
c3 = [0 .4 0];
color_wheel = {c1,c2,c3,'m','b','c','k','y'};
rep_wheel = ceil(n_plots/length(color_wheel));
if rep_wheel>1
    color_wheel = repmat(color_wheel,1,rep_wheel);
end
hold on
use_ste = 1;
for p = 1:n_plots
    data = varargin{p};
    data = squeeze(data);
    thisX = 1:size(data,2);
    thisMu = nanmean(data,1);
    if use_ste
        thisEr = nanstd(data)/sqrt(size(data,1));
    else
        thisEr = nanstd(data);
    end
    thisColor = color_wheel{p};
    
    shadedErrorbar2(thisMu, thisEr,thisX, inputColor,0,alpha);
end
if nargin == 0 %- set to 1 to turn on
   [~,p] = ttest(data);
   if any(p<.05), plot(find(p<.05),0,'ro','MarkerSize',20,'MarkerFaceColor','r'), end
   if any(p<.01), plot(find(p<.01),0,'co','MarkerSize',20,'MarkerFaceColor','c'), end
   if any(p<.005), plot(find(p<.01),0,'go','MarkerSize',20,'MarkerFaceColor','g'), end
   if any(p<.001), plot(find(p<.01),0,'mo','MarkerSize',20,'MarkerFaceColor','m'), end

end

%hold off

function shadedErrorbar2(thisMu, thisEr,thisX, thisColor, HIDE_LINE,alpha)
% function [hLine hShade] = shadedErrorbar2(thisMu, thisEr,thisX, thisColor, HIDE_LINE)
% OG by jw
% made better by ts

if ~exist('alpha','var')
    alpha = .2;
end

if ~exist('thisX','var') || isempty(thisX)
    thisX = 1:length(thisMu);
end
if ~exist('HIDE_LINE','var')
    HIDE_LINE = 0;
end

if ~exist('thisColor','var')
    thisColor = 'b';
end

TURN_OFF_LEGEND = 1; % only one entry per line on legends
if thisMu==0, thisMu = zeros(size(thisX)); end;
thisMu = squeeze(thisMu);
thisEr = squeeze(thisEr);


if size(thisX,1)>size(thisX,2),thisX=thisX'; end
if size(thisMu,1)>size(thisMu,2),thisMu=thisMu'; end
if size(thisEr,1)>size(thisEr,2),thisEr=thisEr'; end
thisNaN = find(~isnan(thisMu.*thisEr));

hLine = plot(thisX, thisMu, 'k-'); hold on
set(hLine,'Color',thisColor,'LineWidth',2);
if HIDE_LINE, set(hLine,'LineStyle','none'); end %- comment this out to show the line
thisX = thisX(thisNaN);
thisMu = thisMu(thisNaN);
thisEr = thisEr(thisNaN);
iFill = [length(thisX):-1:1];
fillX = [thisX thisX(iFill)];
fillY = [thisMu+thisEr thisMu(iFill)-thisEr(iFill)];
hShade = fill(fillX,fillY,thisColor); hold on
set(hShade,'LineStyle','none','FaceAlpha',alpha);     %- alpha values: 0.5
if TURN_OFF_LEGEND
    hAnno = get(hShade,'Annotation');
    hLegEntry = get(hAnno,'LegendInformation');
    set(hLegEntry,'IconDisplayStyle','off')
end
