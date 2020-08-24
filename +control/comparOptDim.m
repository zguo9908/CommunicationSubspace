thetaOpt = [];
deltaOpt = [];
rippleOpt = [];

thetaConOpt = [];
deltaConOpt = [];
rippleConOpt = [];


for j = 1:10
    for i = 2:2
        thetaOpt = [thetaOpt, Patterns(j,i,1).rankRegress.optDimReducedRankRegress];
        deltaOpt = [deltaOpt, Patterns(j,i,2).rankRegress.optDimReducedRankRegress];
        rippleOpt = [rippleOpt, Patterns(j,i,3).rankRegress.optDimReducedRankRegress];
        thetaConOpt = [thetaConOpt, Patterns(j,i,4).rankRegress.optDimReducedRankRegress];
        deltaConOpt = [deltaConOpt, Patterns(j,i,5).rankRegress.optDimReducedRankRegress];
        rippleConOpt = [rippleConOpt, Patterns(j,i,6).rankRegress.optDimReducedRankRegress];

    end
end
result = [mean(thetaOpt), mean(deltaOpt), mean(rippleOpt), mean(thetaConOpt), mean(deltaConOpt), mean(rippleConOpt)]

load("optDims.mat")
resulttable =[ resulttable; table(result)]
save("optDims", "resulttable")