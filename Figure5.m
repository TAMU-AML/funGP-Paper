%% Code for generating Figure 5 in the Appendix

addpath('./ExperimentsSourceFiles/');
rng(1);
[V0, T0] = meshgrid(0.002:0.0002:0.010,340:0.5:360);
Z1 = zeros(size(V0,1),size(V0,2));
M  = 45;
S  = 0.010;
%V0 = V0;
k  = 2000;
P0 = 100000;
Ta = 292;
%T0 = T0;
for i = 1:size(Z1,1)
    for j = 1:size(Z1,2)
        Z1(i,j) = piston([M,S,V0(i,j),k,P0,Ta,T0(i,j)]);
    end
end
surf(V0,T0,Z1)
xlim([0.002,0.01]);
ylim([340,360]);
zlim([0.3,0.7]);
xlabel('V_{0}','FontSize',18);
ylabel('T_{0}','FontSize',18);
zlabel('f(x)','FontSize',18);
set(gca,'FontSize',16)
set(gcf,'PaperSize', [6 5]);
print -painters -dpng Figure5a.png

Z2 = zeros(size(V0,1),size(V0,2));
for i = 1:size(Z2,1)
    for j = 1:size(Z2,2)
        Z2(i,j) = piston([M,S,V0(i,j),k+500,P0,Ta,T0(i,j)]);
    end
end

surf(V0,T0,Z2)
xlim([0.002,0.01]);
ylim([340,360]);
zlim([0.3,0.7]);
xlabel('V_{0}','FontSize',18);
ylabel('T_{0}','FontSize',18);
zlabel('g(x)','FontSize',18);
set(gca,'FontSize',16)
set(gcf,'PaperSize', [6 5]);
print -painters -dpng Figure5b.png

rng(1)
V01 = 0.002 + (0.010 - 0.002).*rand(100,1);
T01 = 340 + (360 - 340).*rand(100,1);

V02 = 0.002 + (0.010 - 0.002).*rand(100,1);
T02 = 340 + (360 - 340).*rand(100,1);

Y1 = zeros(length(V01),1);
for i = 1:length(V01)
    Y1(i) = piston([M,S,V01(i),k,P0,Ta,T01(i)])+ normrnd(0,0.05,1,1);
end

Y2 = zeros(length(V01),1);
for i = 1:length(V01)
    Y2(i) = piston([M,S,V02(i),k+500,P0,Ta,T02(i)])+ normrnd(0,0.05,1,1);
end

hold off
plot(V01,Y1,'o')
hold on
plot(V02,Y2,'x','LineWidth',1.5)
hold off
set(gca,'FontSize',16)
xlabel('V_{0}','FontSize',18);
ylim([0,1]);
%ylabel('y1 and y2','FontSize',18);
legend('y_{1}','y_{2}','Location','northwest')
set(gcf,'PaperSize', [6 5]);
print -painters -dpdf Figure5c.pdf

hold off
plot(T01,Y1,'o')
hold on
plot(T02,Y2,'x','LineWidth',1.5)
hold off
set(gca,'FontSize',16)
xlabel('T_{0}','FontSize',18);
ylim([0,1]);
%ylabel('y1 and y2','FontSize',18);
legend('y_{1}','y_{2}','Location','northwest')
set(gcf,'PaperSize', [6 5]);
print -painters -dpdf Figure5d.pdf
clear