function [ my_str ] = my_sprintf( pattern, varargin )
%MY_SPRINTF Prints nicely and replaces the . in the noise_level by a _

radical = pattern(1:end-4);
extension = pattern(end-3:end);

my_str = [...
    strrep(sprintf(radical, varargin{:}), '.', '_')...
    extension...
    ];

end

