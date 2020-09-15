function Z = offspring(X,s)
 N = numel(X(:,1));                            % population size 
 nMut = sum(X,2);                              % individual mutation load
 sPower = (1-s).^nMut;
 mFitness = mean(sPower);
 fitVec = sPower/mFitness;       % relative fitness\
 y= randsample(N,N,true,fitVec); 
 Z = X(y,:);
end
