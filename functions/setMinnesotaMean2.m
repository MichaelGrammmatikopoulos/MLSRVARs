%% mean for Minnesota prior: zero (diff) or RW (level)

minnesotaPriorMean = NaN(N,1);
for n = 1 : N
    switch tcode(n)
        case {1,4}
            minnesotaPriorMean(n) = 1;
        otherwise
            minnesotaPriorMean(n) = 0;
    end
end
