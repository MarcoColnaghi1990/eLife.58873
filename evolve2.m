%%%. Simulates Muller's ratchet in an asexual population or in a population
%%%. that undergoes LGT, calculating the number of mutations that reach
%%%. fixation.
%%%.
%%%. PARAMETERS:
%%%. nGen = # generations
%%%. U = genome wide mut rate / generation
%%%. s = strength of selection against deleteriouos mutation
%%%. l = LGT rate / individual / generation
%%%. L = eDNA length

function n_fixed_mutations = evolve2(M,nGen,U,s,l,L,additive)

    for t=1:nGen
        M = mutate(M,U,additive);            % introduce new  mutations
        oldMat=M;
        M = offspring(M,s);                  % next generation
    
        if (l>0 && L>0)
            M = LGT(M,l,L,oldMat);           % lateral gene transfer
        end
    
        n_fixed_mutations = sum(min(M));
    end
end


