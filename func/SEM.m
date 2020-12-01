function sem=SEM(X)
    sem=nanstd(X)/sqrt(sum(~isnan(X(:,1)))-1);
end