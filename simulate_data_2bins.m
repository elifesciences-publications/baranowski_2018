
% use this look to simulate data after Bins_analysis_.op
function [stats,Locus_filtered] = simulate_data_2bins(Bins_normbyGene,NumberMutants_exp,NumberMutants_obs,Bins_median,Locus,kmax)

%%this function simulates compares distrubtions of hits in different bins
%%to what might be expected by chance.

BinLoc = [-2.32,2.32]; %this is the location

% take out NaN in Bins_normbyGene
 B = Bins_normbyGene;
 B(any(isnan(B),2),:)=[];

%to make hypothesis testing, use the median and the average number of
%mutants

% initialize matrices
% kmax is number of simulations to run
stats = zeros(length(B),6);
p1 = zeros(1,kmax);
p2 = zeros(1,kmax);
M = zeros(1,kmax);
S = zeros(1,kmax);


% then use the number of mutants for each gene to simulate data call
% compare to simulated wt data.

for i = 1:length(B);
    
    simulated_wt = []; %reset every loop
    simulated_data = []; %reset every loop
    
    for k =1:kmax;
        
        simulated_wt = randn(1,round(mean(NumberMutants_obs))); %simulate wt as normal distribution with mean of 0 and std of 1. Use the mean number of mutants observed.
        %         for j = 1:8
        %             simulated_wt = horzcat(simulated_wt,BinLoc(j)+0.1*randn(1,round(mean(NumberMutants_obs)*Bins_median(1,j))));
        %         end
        
        
        mutant = B(i,1);
        
        %n(i,1) = NumberMutants(mutant);
        n_obs = NumberMutants_obs(mutant);
        n_exp = NumberMutants_exp(mutant);
        
        
        simulated_data = [];
        
        for j = 1:2
            
            simulated_data = horzcat(simulated_data,BinLoc(j)+0.2*randn(1,round(B(i,(j+1))*n_obs)));
            %simulated_data = horzcat(simulated_data,ones(1,round(B(i,(j+1))))*BinLoc(j));
        end
        
        simulated_data1 =simulated_data - median(simulated_data);
        simulated_wt1 = simulated_wt - median(simulated_wt);
        
        if length(simulated_data)>1;
            
            [p1(k) h] = ranksum(simulated_data,simulated_wt);
            
        else
            p1(k) = NaN;
            
            
        end
        
        M(k) = mean(simulated_data);
        
    end
    
    
    
    stats(i,1:6) = horzcat(mutant,n_obs,n_exp,mean(M),mean(p1),std(p1));
    
end

Locus_filtered = Locus(stats(:,1));
%sortedbystd = sortrows(stats,5);
%sortedbyvar = sortrows(stats,8);
%sortedbyttest = sortrows(stats,6);
end
