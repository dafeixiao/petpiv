function write_vel_file(filename, plotline1, plotline2v, plotline3, velData)
    fidv = fopen(filename, 'wt+');
    fprintf(fidv, '%s\n', plotline1);
    fprintf(fidv, '%s\n', plotline2v);
    fprintf(fidv, '%s\n', plotline3);
    fprintf(fidv, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n', velData');
    fclose(fidv);
end
