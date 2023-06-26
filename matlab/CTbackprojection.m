function [ img ] = CTbackprojection( proj, param )
%CTBACKPROJECTION Summary of this function goes here
%   Detailed explanation goes here
disp("hello i am in ct backprojection");
img = zeros(param.nx, param.ny, param.nz, 'single');
for i = 1:param.nProj
    img = img + backprojection(proj(:,:,i),param,i);
    disp(i);
    %tel=img(:,:,2);
    %filename = fullfile('/Users/hi/Downloads/CBCT_Kyungsang_matlab_Feb2015 2/my_progress', sprintf('newinput_%03d.raw', i+1));
    %fid = fopen(filename, 'wb');
    %fwrite(fid, tel, 'float32');
    %fclose(fid);
end
end


