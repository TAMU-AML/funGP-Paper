%% Code for generating Figure 6 in the Appendix

addpath('./ExperimentsSourceFiles/');
rng(1);
[rw, r] = meshgrid(0.05:0.002:0.15,100:1000:50000);
Z = zeros(size(rw,1),size(rw,2));
%rw = rw;
%r  = r;
Tu = 78000;
Hu = 1050;
Tl = 84;
Hl = 760;
L  = 1400;
Kw = 11000;
for i = 1:size(Z,1)
    for j = 1:size(Z,2)
        Z(i,j) = borehole([rw(i,j),r(i,j),Tu,Hu,Tl,Hl,L,Kw]);
    end
end
fig = surf(rw,r,Z);
xlim([0.05,0.15]);
ylim([100,50000]);
zlim([0,200]);
xlabel('r_{w}','FontSize',18);
ylabel('r','FontSize',18);
zlabel('f(x)','FontSize',18);
set(gca,'FontSize',18)
set(gcf,'PaperSize', [6 5]);
print -painters -dpng Figure6a.png

for i = 1:size(Z,1)
    for j = 1:size(Z,2)
        Z(i,j) = borehole([rw(i,j),r(i,j),Tu,Hu,Tl,Hl,L+50,Kw]);
    end
end
fig = surf(rw,r,Z);
xlim([0.05,0.15]);
ylim([100,50000]);
zlim([0,200]);
xlabel('r_{w}','FontSize',18);
ylabel('r','FontSize',18);
zlabel('g(x)','FontSize',18);
set(gca,'FontSize',18)
set(gcf,'PaperSize', [6 5]);
print -painters -dpng Figure6b.png



rw01 = 0.05 + (0.15 - 0.05).*rand(100,1);
r01  = 100 + (50000 - 100).*rand(100,1);
Tu = 78000;
Hu = 1050;
Tl = 84;
Hl = 760;
L  = 1400;
Kw = 11000;
x1 = [rw01,r01];
Y1 = zeros(length(rw01),1);
for i = 1:length(rw01)
    Y1(i) = borehole([rw01(i),r01(i),Tu,Hu,Tl,Hl,L,Kw])+ normrnd(0,10,1,1);
end

rw02 = 0.05 + (0.15 - 0.05).*rand(100,1);
r02  = 100 + (50000 - 100).*rand(100,1);
x2 = [rw02,r02];
Y2 = zeros(length(rw01),1);
for i = 1:length(rw01)
    Y2(i) = borehole([rw02(i),r02(i),Tu,Hu,Tl,Hl,L+50,Kw])+ normrnd(0,10,1,1);
end

fig = plot(rw01,Y1,'o');
hold on
plot(rw02,Y2,'x','LineWidth',1.5)
hold off
set(gca,'FontSize',18)
xlabel('r_{w}','FontSize',18);
ylim([0,200]);
%ylabel('y1 and y2','FontSize',18);
legend('y_{1}','y_{2}','Location','northwest')
set(gcf,'PaperSize', [6 5]);
print -painters -dpdf Figure6c.pdf


fig = plot(r01,Y1,'o');
hold on
plot(r02,Y2,'x','LineWidth',1.5)
hold off
set(gca,'FontSize',18)
xlabel('r','FontSize',18);
ylim([0,200]);
%ylabel('y1 and y2','FontSize',18);
legend('y_{1}','y_{2}','Location','northwest')
set(gcf,'PaperSize', [6 5]);
print -painters -dpdf Figure6d.pdf
clear