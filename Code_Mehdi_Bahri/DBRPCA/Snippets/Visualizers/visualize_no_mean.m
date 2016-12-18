function [ ] = visualize_no_mean( vars, params )
%VISUALIZE_NO_MEAN Summary of this function goes here
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

subplot(2,4,1), imshow(vars.X(:,:,1), [])
subplot(2,4,2), imshow(Uc*T(:,:,1)*Ur', [])
subplot(2,4,3), imshow(vars.E(:,:,1), [])

% subplot(1,5,1), imshow(X(:,:,1))
% subplot(1,5,2), imshow(vars.Uc*vars.T(:,:,1)*vars.Ur')
% subplot(1,5,3), imshow(vars.E(:,:,1))

subplot(2,4,4), stemplot(Uc,1,1);
subplot(2,4,5), stemplot(Ur,1,1);
subplot(2,4,6), meshplot(T(:,:,1),2,1);
subplot(2,4,7), stemplot(T(:,:,1),2,1);
drawnow;
    
end

