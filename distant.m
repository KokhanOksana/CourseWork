function d = distant(x, y)
    d = 0;
    for k = 1:size(x)
        d = d + (x(k) - y(k))^2;
    end;
    d = sqrt(d);
end

