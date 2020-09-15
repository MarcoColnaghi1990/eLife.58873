%%

figure;hold on
colors=get(gca,'colororder');

%%%. System Parameters
nGen=10000;                         % Number of generations
N=5000;                             % Population size
g_vector= round(10.^(2:1/6:3));     % Number of genes; 1 gene =~ 1000bp
nRuns = 50;                         % # of Iterations
Lvector = [0 1 5 10 100000];        % eDNA length
l = 0.01;                           % LGT rate / individual / generation
u_bp= 10*10^(-8);                   % Mutation rate per bp 
u = 1000*u_bp;                      % Mutation rate per gene
s = 0.001;                          % Strength of selection agains deleterious mutations
inLoad = 0.0;                       % Initial mutation load
dmdt={[];[];[];[];[];};

%%%. Simulation
for L1 = 1:numel(Lvector)           
    L=Lvector(L1);
    n_fixed_mutations=zeros(numel(g_vector),nRuns);
    
    for g1 = 1:numel(g_vector)
        g = g_vector(g1);
        U=u*g;                                        % Genome-wide mutation rate
        if L1 ==numel(Lvector)
            L=round(g/5);
        end
        disp(['L = ', num2str(L), '. g = ', num2str(g), '     n_0 = ',num2str(1/exp(U/s) *N)])
        
        for r1 = 1:nRuns
            if inLoad>0
                rand_mat = rand(N,g);                 % generate random matrix (for mutations)
                M = double(rand_mat < inLoad);        % N x g matrix; 0 = wild-type, 1 = mutant
            else
                M = zeros(N,g);
            end
            
            n_fixed_mutations(g1,r1) = evolve2(M,nGen,U,s,l,L,true);
            disp(['     run ',num2str(r1),' of ', num2str(nRuns), '     *** ', num2str(n_fixed_mutations(g1,r1)), ' fixed mutations ***'])
        end
    end
    
    n_fixed_mutations = n_fixed_mutations/nGen;
    dmdt{L1} = n_fixed_mutations;
    er = sqrt(var(n_fixed_mutations,0,2));
    y4   = mean(n_fixed_mutations,2);
    
    errorbar(g_vector,y4,er,'linewidth',2,'linestyle',':','marker','o')
    set(gca,'Xscale','log')
    drawnow
end
