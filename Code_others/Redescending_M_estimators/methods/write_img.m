function write_img(T,name)

        gp = im2uint8(T);
    imwrite(gp,name);

return;