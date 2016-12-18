function [ ps ] = extract_psnr( results )
%EXTRACT_PSNR Extracts the list of PSNRs from the measure structs
%
% Mehdi Bahri - Imperial College London
% July, 2016

ps = zeros(1, length(results));

for i = 1:length(results)
    if ~isa(results{i}, 'MException')
        ps(i) = results{i}.MSIQA.psnr;
    else
        ps(i) = 0;
    end
end

end

