function p = get_measure_err(a, measure)

if isa(a, 'MException')
    if strcmp(measure, 'msam') || strcmp(measure, 'rel_norm');
        p = 1;
    else
        p = 0;
    end
else
    p = a.MSIQA_E.(measure);
end

end