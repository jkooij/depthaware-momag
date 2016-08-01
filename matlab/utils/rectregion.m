function r = rectregion(x0, y0, width, height)
    x1 = x0 + width;
    y1 = y0 + height;

    r = struct;
    r.x0 = x0;
    r.y0 = y0;
    r.x1 = x1;
    r.y1 = y1;
    r.width = width;
    r.height = height;
    r.x0y0wh = [x0 y0 width height];
    r.x0x1y0y1 = [x0 x1 y0 y1];
end