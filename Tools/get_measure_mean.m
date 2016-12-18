function p = get_measure_mean(a, measure)

if isa(a, 'MException')
    if strcmp(measure, 'msam') || strcmp(measure, 'rel_norm');
        p = 1;
    else
        p = 0;
    end
else
    p = a.MSIQA_M.(measure);
end

end