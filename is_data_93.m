function [data]=is_data_93(data) 
%% If try to use the data of 93, use this part to correct the data %%          
for a = 1:4
    for b = 1:14
        for c = 1:2
            for d = 1:131072
                if (data(d,c,b,a)) < 0
                    data(d,c,b,a) = data(d,c,b,a) + 256;
                end
            end
        end
    end
end
