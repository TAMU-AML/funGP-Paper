addpath('./ExperimentsSourceFiles/','./Algorithm');
%% Code for generating Figure 1
ntrain = 100;
sigma_f = 5;
sigma_n = 0.5;
theta = 0.2;
rng(1);
x = linspace(0,1,1000)';
CovMat_F = ComputeCovMat(x,x,theta,"SqExp");
CovMat_F = (sigma_f^2)*CovMat_F;
mu = zeros(1000,1);
F = mvnrnd(mu,CovMat_F,1)' ;
FSq = F.^2;
l2norm = sqrt(trapz(x,FSq));
perturbationRegion = [0.2,0.8];
fun = @(x) 0.33*sin(pi*(x-perturbationRegion(1))/(perturbationRegion(2)-perturbationRegion(1)));
perturb_index = find(x>=perturbationRegion(1) & x<=perturbationRegion(2));
Gx = F;
Gx(perturb_index) = Gx(perturb_index)+ fun(x(perturb_index));
hold off;
plot(x,F,'LineWidth',1.5);
hold on;
plot(x,Gx,'--','LineWidth',1.5);
xlabel('x');
legend('f(x)','g(x)','Location','northwest')
set(gca,'FontSize',18)
set(gcf,'PaperSize', [6 5]);
hold off;
print -painters -dpdf Figure1a.pdf

index1 = randsample(1:1000,ntrain)';
index2 = randsample(1:1000,ntrain)';
y1 = F(index1) + normrnd(0,sigma_n,ntrain,1);
y2 = Gx(index2) + normrnd(0,sigma_n,ntrain,1);
p = [];
hold off;
p(1) = plot(x(index1),y1,'o');
hold on;
p(2) = plot(x(index2),y2,'x');
xlabel('x');
legend([p(1),p(2)],'y_{1}','y_{2}','Location','northwest')
set(gca,'FontSize',18)
set(gcf,'PaperSize', [6 5]);
hold off;
print -painters -dpdf Figure1b.pdf
clear



function[CovMat] = ComputeCovMat(x1,x2,theta,CovType)
    mat_size1 = size(x1,1);
    mat_size2 = size(x2,1);
    CovMat = zeros(mat_size1,mat_size2);
    if CovType == "SqExp"
        for i = 1:mat_size1
                for j = 1:mat_size2
                   t = 0;
                    for p = 1:length(theta)
                      t =  t + ((x1(i,p) - x2(j,p))/theta(p))^2;
                    end
                    CovMat(i,j) =   exp(-0.5*t);
                end
        end
    end
    
    if CovType == "Matern"
        for i = 1:mat_size1
                for j = 1:mat_size2
                   t = 0;
                    for p = 1:length(theta)
                      t =  t + ((x1(i,p) - x2(j,p))/theta(p))^2;
                    end
                    t = sqrt(t);
                    CovMat(i,j) =   (1+(sqrt(3)*t))*exp(-sqrt(3)*t);
                end
        end
    end
    
end