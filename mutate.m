
%%% Mutation function. It draws a random integer X from a Poisson
%%% distribution for each individual, which is the number of new mutations
%%% acquired. For each individual, it selects X loci in its genome which,
%%% if unmutated, acquire a new mutation.
%%% If additive == True, the mutations add up, otherwise only 1
%%% mutations/locus

function X = mutate(X,U,additive)
 N = numel(X(:,1));                            % population size 
 g = numel(X(1,:));                            % genome size 
 mutString = random('Poisson',U,N,1);

 nMut = sum(mutString);
 if nMut>0
    yStr = randi([1,g],1,nMut); a=1;
    pStr = zeros (1,nMut);
    mut=find(mutString>0);
    for k=1:numel(mut)
    	i=mut(k);
    	x=repmat(i,1,mutString(i));
    	pStr(a:a+mutString(i)-1)=x;
        a=a+mutString(i);
    end
    z= sum ([g.*(pStr-1);yStr]);
    X(z) = X(z)+1;
 end
 
 if ~additive
    X = min(X,1);
 end
 
end
