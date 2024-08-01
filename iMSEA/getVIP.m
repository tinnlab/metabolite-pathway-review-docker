function vipScore = getVIP(X,y,ncomp)
%GETVIP 此处显示有关此函数的摘要
%   此处显示详细说明
% X = normalize(X);
X = normalize(X,'center','mean');
y = normalize(y);

[XL,yl,XS,~,~,~,~,stats] = plsregress(X,y,ncomp);
W0 = stats.W ./ sqrt(sum(stats.W.^2, 1));% Calculate the normalized PLS weights
p = size(XL, 1);
sumSq = sum(XS.^2, 1).*sum(yl.^2, 1); % 按行求和后点乘(计算相关系数)
vipScore = sqrt(p * sum(sumSq.*(W0.^2), 2) ./ sum(sumSq, 2));
% indVIP = find(vipScore >= 1);
% scatter(1:length(vipScore),vipScore,'x')
% hold on
% scatter(indVIP,vipScore(indVIP),'rx')
% plot([1 length(vipScore)],[1 1],'--k')
% hold off
% axis tight
% xlabel('Predictor Variables')
% ylabel('VIP Scores')
end

