close all; clear; clc;

y = 5;
x = 5;
xLimit = [0, 6];
lineLength = 0.08;
ylim([y-1, y+1])
xlim(xLimit)
line(xLimit, [y,y], "LineWidth", 2)
line([1,1], [y-lineLength, y+lineLength], "LineWidth", 2)

txt = {'T_{m-2}','T_{m-1}','T_m','T_{m+1}','T_{m+2}'};

for i = 1:5
    line([i,i], [y-lineLength, y+lineLength], "LineWidth", 2)
    text(i-0.05, y-0.2, txt(i))
end

TM = 3.7;
line([TM,TM], [y-lineLength, y+lineLength], "LineWidth", 2)
    text(TM-0.05, y+0.18, "T_M")

set(gca,'YTick',[]);
set(gca,'XTick',[]);

x = [1,2,3];
y = [5.8,5.6,5.2];
p = polyfit(x,y,2);

x1 = xLimit(1):0.1:xLimit(2);
y1 = polyval(p,x1);

hold on
plot(x1, y1, "LineWidth", 2)
line([TM-0.5,TM+0.5], [polyval(p,TM)+0.32, polyval(p,TM)-0.32], "LineWidth", 2, 'Color', 'red')


h1 = plot(x, y, 'o', 'Color', 'black');
set(h1, 'markerfacecolor', get(h1, 'color'));
h2 = plot(TM, polyval(p,TM), 'o', 'Color', 'red');
set(h2, 'markerfacecolor', get(h2, 'color'));


% plot 2 (warning: lots of duplicate code...)
figure
y = 5;
x = 5;
xLimit = [0, 6];
lineLength = 0.08;
ylim([y-1, y+1])
xlim(xLimit)
line(xLimit, [y,y], "LineWidth", 2)
% line([1,1], [y-lineLength, y+lineLength], "LineWidth", 2)

txt = {'T_{m-2}','T_{m-1}','T_m','T_{m+1}','T_{m+2}'};

for i = 1:5
    line([i,i], [y-lineLength, y+lineLength], "LineWidth", 2)
    text(i-0.05, y-0.2, txt(i))
end

TM = 4.3;
line([TM,TM], [y-lineLength, y+lineLength], "LineWidth", 2)
    text(TM-0.05, y+0.18, "T_M")

set(gca,'YTick',[]);
set(gca,'XTick',[]);

x = [2,3,TM];
y = [5.8,5.5,4.8];
p = polyfit(x,y,2);

x1 = xLimit(1):0.1:xLimit(2);
y1 = polyval(p,x1);

x2 = polyder(p);
y2 = polyval(p,x2);

hold on
plot(x1, y1, "LineWidth", 2)
% line([4-0.5,4+0.5], [5.3, 4.7], "LineWidth", 2, 'Color', 'red')


h1 = plot(x, y, 'o', 'Color', 'black');
set(h1, 'markerfacecolor', get(h1, 'color'));
h2 = plot(4, polyval(p,4), 'o', 'Color', 'red');
set(h2, 'markerfacecolor', get(h2, 'color'));


a=1;
b=10;
xy = a:b;
[X, Y] = meshgrid(xy);

water = 10;
rx = (b - a) .* rand(b*water, 1) + a;
ry = (b - a) .* rand(b*water, 1) + a;
ry = ry - b;

figure
h1 = plot(X, Y, 'o', 'Color', 'blue', "LineWidth", 0.7);
set(h1, 'markerfacecolor', 'blue');
hold on
h1 = plot(rx, ry,'o','Color', 'blue');
set(h1, 'markerfacecolor', 'blue');


xlim([0,11])
ylim([-11,11])
set(gca,'YTick',[]);
set(gca,'XTick',[]);
axis square


