function [svd_v, svd_a] = run_svd_on_grid( ...
    Xg, Yg, Zg, Region, scale, ...
    Xp, Yp, Zp, Up, Vp, Wp, axp, ayp, azp, T, ...
    delta, nu, scale_unit, Nm_v, Nm_a, weight, UNCUp, UNCVp, UNCWp,cutoff)

    N1 = length(Xg);
    N2 = length(Yg);
    N3 = length(Zg);
    totalN = N1 * N2 * N3;

    svd_v = zeros(totalN, 18); % (x,y,z,u,v,w,dudx,...,dwdt)
    svd_a = zeros(totalN, 18); % (x,y,z,ax,ay,az,...)

% test code, dafei
%     K = zeros(totalN,1);
%     for idx = 1:totalN
%         [~,~,k] = index_to_ijk(idx,N1,N2,N3);
%         K(idx) = k;
%     end
%     
%     unique(K(end-10:end))



    parfor indexcount = 1:totalN
    % for indexcount = 1:totalN
        [i, j, k] = index_to_ijk(indexcount, N1, N2, N3);

        point = [Xg(i), Yg(j), Zg(k)];

        % Velocity SVD (second order)
        datarow_v = zeros(1,18);
        datarow_v(1:3) = point;

        [Sol_v, count_v] = solver2_second(point, Xp, Yp, Zp, Up, Vp, Wp, Region*scale*3, 11, T, UNCUp, UNCVp, UNCWp,cutoff,Zg);
        scale_2 = 1;
        while count_v < 11
            scale_2 = scale_2 + 0.1;
            [Sol_v, count_v] = solver2_second(point, Xp, Yp, Zp, Up, Vp, Wp, Region*scale*3*scale_2, 11, T, UNCUp, UNCVp, UNCWp,cutoff,Zg);
        end

        datarow_v(4:18) = [Sol_v(1:12); Sol_v(end-2:end)];
        svd_v(indexcount,:) = datarow_v;

        % Acceleration SVD (first order in time)
        datarow_a = zeros(1,18);
        datarow_a(1:3) = point;

        [Sol_a, count_a] = solver2_first_time(point, Xp, Yp, Zp, axp, ayp, azp, Region*scale*3, 5, T);
        scale_3 = 1;
        while count_a < 5
            scale_3 = scale_3 + 0.1;
            [Sol_a, count_a] = solver2_first_time(point, Xp, Yp, Zp, axp, ayp, azp, Region*scale*3*scale_3, 5, T);
        end

        if count_a >= 5
            datarow_a(4:end) = Sol_a;
        end
        svd_a(indexcount,:) = datarow_a;
    end
end
