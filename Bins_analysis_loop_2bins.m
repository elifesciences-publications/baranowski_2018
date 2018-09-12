function [Bins_normbyReadCount,Bins_normbyGene,NumberMutants_obs,NumberMutants_exp,Bins_median] = Bins_analysis_loop_2bins(Low,High,control_msm)

% import Bin TA counts and Locus indentifier as cell array
% column 1 is # of hit TAs
% column 2 is # of Total TAs
% column 3 is Fraction TAs hit
% column 4 is Total Reads
% column 5 is Average Reads/TA

%clearvars -except Low1 High1 Low2 High2 Locus control_msm

Bin_cell = {Low, High};
Bins_combined=[[1:9984]']; %makes an array from 1 to 9984
BinSize = [1,1]; %this is the expected bin width
Correction = [1,1];
BinLoc = [-2.32,2.32]; % this is the expected position of each bin as computed for 10% of a normal distribution


%% combine each the Total Reads for each Bin into a large arrary

for i = 1:length(Bin_cell)
    
    Bin = Bin_cell{i};
    Bins_combined = horzcat(Bins_combined,Bin(:,4));
    
end

Bins_normbyReadCount = [Bins_combined(:,1)];

for i = 1:length(Bin_cell)
    
    Bin = Bin_cell{i};
    Bins_normbyReadCount = horzcat(Bins_normbyReadCount,Bins_combined(:,(i+1))/sum(Bin(:,4)));
    
end

%% the number of cells for each mutant assuming you sorted a total of 300K cells/bin

total_reads = sum(Bins_normbyReadCount(:,2:3),2);
NumberMutants_obs = total_reads./sum(total_reads)*300000*2;
total_reads = sum(control_msm);
NumberMutants_exp= control_msm/total_reads*3000002*2;
% NumberMutants=geomean([NumberMutants1,NumberMutants2],2);

%% normalizes each Bin by width of Bin to account for nonlinear Bin size and
% corrects for the unevenness of the bins.

Bins_normbyBinSize = [Bins_combined(:,1)];

for i = 1:length(Bin_cell)
    
    Bin = Bin_cell{i};
    Bins_normbyBinSize = horzcat(Bins_normbyBinSize,Bins_normbyReadCount(:,(i+1))/(BinSize(i)*Correction(i)));
    
end

total = sum(Bins_normbyBinSize(:,2:3),2);
Bins_normbyGene= [Bins_combined(:,1)];

for i = 1:length(Bin_cell)
    
    Bin = Bin_cell{i};
    Bins_normbyGene= horzcat(Bins_normbyGene,Bins_normbyBinSize(:,(i+1))./total);
    
end

Bins_median = nanmedian(Bins_normbyGene(:,2:3),1);

end
