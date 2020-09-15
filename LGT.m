
%%% LGT function. It takes as input the probability of lgt l.
%%% string length L and the old matrix.
%%%
%%% For each individual, generate a random number; if it is smaller than l,
%%% select a random locus x between 1 and g-L+1. Select a random row from the
%%% old matrix. Elements x to x+L-1 of the new matrix become equal to those
%%% of the previous matrix. 

function X = LGT(X,l,L,oldMat)
    N = numel(X(:,1));                            % population size 
    g = numel(X(1,:));                            % genome size 
 
    nLGT = random('Bino',N,l);                    % # of LGT transfer events in the population
    whoLGT = randsample(1:N,nLGT);                 % who undergoes LGT

    if sum(nLGT)>0
        LGTloci = randi([1 g-L+1],1,sum(nLGT));   % locus where rec begins
        if LGTloci+L>g
            diff = LGTloci+L -g-1;
            X(whoLGT, LGTloci:end) = oldMat(whoLGT, LGTloci:end);  
            X(whoLGT, 1:diff) = oldMat(whoLGT, LGTloci:1:diff);    
        else
            X(whoLGT, LGTloci:LGTloci+L-1) = oldMat(whoLGT, LGTloci:LGTloci+L-1);
                                                % recombination 
        end
    end
end