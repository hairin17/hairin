function [corr_val] = vec_corr(vec_j,vec_k)

cor = zeros(2,2)
for i = 1:95
cor = corrcoef(vec_j(i,:),vec_k(i,:))

corr_val(i) = cor(1,2)
corr_val(i) = log((1+corr_val(i))./(1-corr_val(i)))/2;

end
% corr_val(i) = corr_val'
end

