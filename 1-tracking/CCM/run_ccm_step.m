function [ccm_v, ccm_a] = run_ccm_step( ...
    Xg, Yg, Zg, Region, scale, ...
    Xp, Yp, Zp, Up, Vp, Wp, axp, ayp, azp, T, ...
    UNCUp, UNCVp, UNCWp, UNCAxp, UNCAyp, UNCAzp, ...
    svd_v, svd_a, ...
    delta, nu, scale_unit, ...
    Nm_v, Nm_a, weight)

    N1 = length(Xg);
    N2 = length(Yg);
    N3 = length(Zg);
    totalN = N1 * N2 * N3;

    ccm_v = svd_v;
    ccm_a = svd_a;

    % Reshape SVD results into 3D fields
    uprior  = reshape(svd_v(:,4), N1, N2, N3);
    vprior  = reshape(svd_v(:,5), N1, N2, N3);
    wprior  = reshape(svd_v(:,6), N1, N2, N3);
    uprior2 = reshape(svd_a(:,4), N1, N2, N3);
    vprior2 = reshape(svd_a(:,5), N1, N2, N3);
    wprior2 = reshape(svd_a(:,6), N1, N2, N3);

    dudx  = reshape(svd_v(:,7),  N1, N2, N3);
    dudy  = reshape(svd_v(:,8),  N1, N2, N3);
    dudz  = reshape(svd_v(:,9),  N1, N2, N3);
    dvdx  = reshape(svd_v(:,10), N1, N2, N3);
    dvdy  = reshape(svd_v(:,11), N1, N2, N3);
    dvdz  = reshape(svd_v(:,12), N1, N2, N3);
    dwdx  = reshape(svd_v(:,13), N1, N2, N3);
    dwdy  = reshape(svd_v(:,14), N1, N2, N3);
    dwdz  = reshape(svd_v(:,15), N1, N2, N3);
    dudt  = reshape(svd_v(:,16), N1, N2, N3);
    dvdt  = reshape(svd_v(:,17), N1, N2, N3);
    dwdt  = reshape(svd_v(:,18), N1, N2, N3);

    dudx2 = reshape(svd_a(:,7),  N1, N2, N3);
    dudy2 = reshape(svd_a(:,8),  N1, N2, N3);
    dudz2 = reshape(svd_a(:,9),  N1, N2, N3);
    dvdx2 = reshape(svd_a(:,10), N1, N2, N3);
    dvdy2 = reshape(svd_a(:,11), N1, N2, N3);
    dvdz2 = reshape(svd_a(:,12), N1, N2, N3);
    dwdx2 = reshape(svd_a(:,13), N1, N2, N3);
    dwdy2 = reshape(svd_a(:,14), N1, N2, N3);
    dwdz2 = reshape(svd_a(:,15), N1, N2, N3);
    dudt2 = reshape(svd_a(:,16), N1, N2, N3);
    dvdt2 = reshape(svd_a(:,17), N1, N2, N3);
    dwdt2 = reshape(svd_a(:,18), N1, N2, N3);

    % Gradients for viscous terms (non-uniform grid)
    [dudxdy, dudxdx, dudxdz] = grad3nu(dudx, Xg(:), Yg(:), Zg(:));
    [dudydy, dudydx, dudydz] = grad3nu(dudy, Xg(:), Yg(:), Zg(:));
    [dudzdy, dudzdx, dudzdz] = grad3nu(dudz, Xg(:), Yg(:), Zg(:));

    [dvdxdy, dvdxdx, dvdxdz] = grad3nu(dvdx, Xg(:), Yg(:), Zg(:));
    [dvdydy, dvdydx, dvdydz] = grad3nu(dvdy, Xg(:), Yg(:), Zg(:));
    [dvdzdy, dvdzdx, dvdzdz] = grad3nu(dvdz, Xg(:), Yg(:), Zg(:));

    [dwdxdy, dwdxdx, dwdxdz] = grad3nu(dwdx, Xg(:), Yg(:), Zg(:));
    [dwdydy, dwdydx, dwdydz] = grad3nu(dwdy, Xg(:), Yg(:), Zg(:));
    [dwdzdy, dwdzdx, dwdzdz] = grad3nu(dwdz, Xg(:), Yg(:), Zg(:));

    visc_d = zeros(totalN, 3);
    visc_d(:,1) = reshape(-nu*(dudxdx + dudydy + dudzdz), N1*N2*N3, 1);
    visc_d(:,2) = reshape(-nu*(dvdxdx + dvdydy + dvdzdz), N1*N2*N3, 1);
    visc_d(:,3) = reshape(-nu*(dwdxdx + dwdydy + dwdzdz), N1*N2*N3, 1);

    % Curl of viscous terms
    [dudyv, dudxv, dudzv] = grad3nu(-nu*(dudxdx + dudydy + dudzdz), Xg(:), Yg(:), Zg(:));
    [dvdyv, dvdxv, dvdzv] = grad3nu(-nu*(dvdxdx + dvdydy + dvdzdz), Xg(:), Yg(:), Zg(:));
    [dwdyv, dwdxv, dwdzv] = grad3nu(-nu*(dwdxdx + dwdydy + dwdzdz), Xg(:), Yg(:), Zg(:));

    curl = zeros(totalN,3);
    curl(:,1) = -reshape(dwdyv - dvdzv, N1*N2*N3, 1);
    curl(:,2) = -reshape(dudzv - dwdxv, N1*N2*N3, 1);
    curl(:,3) = -reshape(dvdxv - dudyv, N1*N2*N3, 1);
    curl = curl ./ (scale_unit^2);

    % Main CCM loop over grid points
    parfor indexcount = 1:totalN
        [i, j, k] = index_to_ijk(indexcount, N1, N2, N3);
        point = [Xg(i), Yg(j), Zg(k)];

        % Local neighborhood indices
        [is, ie, js, je, ks, ke] = local_neighborhood(i, j, k, N1, N2, N3, 3);

        % Build S and S2 matrices (covariance data)
        [s_mat, s2_mat] = build_S_matrices( ...
            uprior, uprior2, ...
            vprior, vprior2, ...
            wprior, wprior2, ...
            dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, ...
            dudx2, dudy2, dudz2, dvdx2, dvdy2, dvdz2, dwdx2, dwdy2, dwdz2, ...
            dudt, dvdt, dwdt, ...
            dudt2, dvdt2, dwdt2, ...
            is, ie, js, je, ks, ke);

        % Inverse covariance (or diagonal fallback)
        s  = build_inverse_covariance(s_mat, 33);
        s2 = build_inverse_covariance(s2_mat, 15);

        % Prior state for velocity
        xprior = zeros(33,1);
        xprior(1:3)   = [svd_v(indexcount,4); svd_v(indexcount,5); svd_v(indexcount,6)];
        xprior(4:12)  = [dudx(i,j,k); dudy(i,j,k); dudz(i,j,k); ...
                         dvdx(i,j,k); dvdy(i,j,k); dvdz(i,j,k); ...
                         dwdx(i,j,k); dwdy(i,j,k); dwdz(i,j,k)];
        xprior(13:30) = [dudxdx(i,j,k), dudydy(i,j,k), dudzdz(i,j,k), ...
                         dudxdy(i,j,k), dudxdz(i,j,k), dudydz(i,j,k), ...
                         dvdxdx(i,j,k), dvdydy(i,j,k), dvdzdz(i,j,k), ...
                         dvdxdy(i,j,k), dvdxdz(i,j,k), dvdydz(i,j,k), ...
                         dwdxdx(i,j,k), dwdydy(i,j,k), dwdzdz(i,j,k), ...
                         dwdxdy(i,j,k), dwdxdz(i,j,k), dwdydz(i,j,k)];
        xprior(31:33) = [dudt(i,j,k); dvdt(i,j,k); dwdt(i,j,k)];

        up = zeros(15,1);

        % Velocity CCM
        [Sol_v, count_v] = solver3_2(point, Xp, Yp, Zp, Up, Vp, Wp, ...
                                     UNCUp, UNCVp, UNCWp, ...
                                     Region*scale*2, Nm_v, xprior, s, up, ...
                                     weight, T, 0.01, 4);

        scale_2 = 1;
        while count_v < Nm_v
            scale_2 = scale_2 + 0.1;
            [Sol_v, count_v] = solver3_2(point, Xp, Yp, Zp, Up, Vp, Wp, ...
                                         UNCUp, UNCVp, UNCWp, ...
                                         Region*scale*2*scale_2, Nm_v, xprior, s, up, ...
                                         weight, T, 0.01, 4);
        end

        datarow_v = zeros(1,18);
        datarow_v(1:3) = point;
        if count_v >= 1
            datarow_v(1:18) = [point, Sol_v(1:15)'];
        end
        datarow_v(7:15) = datarow_v(7:15) ./ scale_unit;
        ccm_v(indexcount,:) = datarow_v;

        % Acceleration CCM
        xprior2 = zeros(15,1);
        xprior2(1:3)  = [uprior2(i,j,k); vprior2(i,j,k); wprior2(i,j,k)];
        xprior2(4:12) = [dudx2(i,j,k); dudy2(i,j,k); dudz2(i,j,k); ...
                         dvdx2(i,j,k); dvdy2(i,j,k); dvdz2(i,j,k); ...
                         dwdx2(i,j,k); dwdy2(i,j,k); dwdz2(i,j,k)];
        xprior2(13:15) = [dudt2(i,j,k); dvdt2(i,j,k); dwdt2(i,j,k)];

        up = zeros(15,1);

        [Sol_a, count_a] = solver3_acc_2(point, Xp, Yp, Zp, axp, ayp, azp, ...
                                         UNCAxp, UNCAyp, UNCAzp, ...
                                         Region*2*scale, Nm_a, xprior2, s2, up, ...
                                         weight, curl(indexcount,:), T, 0.01, 4);

        scale_3 = 1;
        while count_a < Nm_a
            scale_3 = scale_3 + 0.1;
            [Sol_a, count_a] = solver3_acc_2(point, Xp, Yp, Zp, axp, ayp, azp, ...
                                             UNCAxp, UNCAyp, UNCAzp, ...
                                             Region*2*scale*scale_3, Nm_a, ...
                                             xprior2, s2, up, weight, ...
                                             curl(indexcount,:), T, 0.01, 4);
        end

        datarow_a = zeros(1,18);
        datarow_a(1:3) = point;
        if count_a >= 1
            datarow_a(1:18) = [point, Sol_a'];
        end
        ccm_a(indexcount,:) = datarow_a;
    end
end
