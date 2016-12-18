function [ ] = visualize_with_mean( vars, params )
%VIZUALIZE_WITH_MEAN Summary of this function goes here
%   Detailed explanation goes here

if isfield(vars, 'A')
    Uc = vars.A;
else
    Uc = vars.Uc;
end

if isfield(vars, 'B')
    Ur = vars.B;
else
    Ur = vars.Ur;
end

if isfield(vars, 'K')
    T = vars.K;
else
    T = vars.T;
end

subplot(2,4,1); imshow(Uc*T(:,:,1)*Ur', []); title('Rec. L');
subplot(2,4,2); imshow(vars.M(:,:,1), []); title('Rec. M');
subplot(2,4,3); imshow(Uc*T(:,:,1)*Ur' + vars.M(:,:,1), []); title('Rec. O');
subplot(2,4,4); imshow(Uc*T(:,:,1)*Ur' + vars.M(:,:,1) + vars.E(:,:,1), []); title('Rec. X');
subplot(2,4,5); imshow(vars.E(:,:,1), []); title('Rec. E');
subplot(2,4,6); stemplot(Uc,1,1);
subplot(2,4,7); stemplot(Ur,1,1);
if isfield(params, 'ground_O')
    subplot(2,4,8); imshow(params.g_O(:,:,1), []); title('Original')
end
drawnow;
    
end

