%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DESCRIPTION:
%   script for funGP hypothesis test for wind turbine datasets. Must be called along with definition for turbine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(1);

if turbine=="WT3" || turbine=="WT4"
 type = "Offshore";
 years = ["2007","2008","2009","2010"];
else
 type = "Inland";
 years = ["2008","2009","2010","2011"];
end

nYears = length(years);
yearlyData = cell(nYears,1);
wscol = 1; 
ycol = 7;
confLevel = 0.9;
nCases = nchoosek(nYears,2);
percentPoints = zeros(1,nCases);
xtest = linspace(3,15,1000)';
muDiff = zeros(size(xtest,1),nCases);
band = zeros(size(xtest,1),nCases);

for i = 1:nYears
    yearlyData{i} = dlmread(strcat(type," Wind Farm Dataset2(",turbine,")_match_",years(i),".txt"),' ',1,2);
end
idx = 0;
for i = 1:(nYears-1)
    for j = (i+1):nYears
        idx = idx + 1;
        x1 = yearlyData{i}(:,wscol);
        y1 = yearlyData{i}(:,ycol);
        x2 = yearlyData{j}(:,wscol);
        y2 = yearlyData{j}(:,ycol);
        result = funGP(x1,y1,x2,y2,xtest,confLevel);
        percentPoints(idx) = result.nPoints*100/size(xtest,1);
        muDiff(:,idx) = result.muDiff;
        band(:,idx) = result.band(:,1);
    end
end
fprintf('Percentage of statistically different test points for %s:\n',turbine);
disp(percentPoints);




