function ensure_dir(d)
    if ~exist(d, 'dir')
        mkdir(d);
    end
end
