% generated a random tensor that has zeros at some index of leve missing_ratio, while 1 at other
% index 
function f = get_harda_oper( size_tens , missing_ratio )
%GET_HARDA_OPER Summary of this function goes here
%   Detailed explanation goes here

    ss = size(size_tens);
    ss = ss(2);
    
    n= 1;
    
    for i = 1: ss
        n = n * size_tens(i);
    end
 
    tmp = randperm( n ) ;
    f = zeros ( n , 1 ) ;
    
    s = (1-missing_ratio)*n;
    s = floor(s) + 1;
    if s > n
        s = n;
    end
    
    f ( tmp ( 1 : s ) ) = ones( s , 1 ) ;
 
    f = vec2tens(f,size_tens,1:ss);

end

