function [s_mat, s2_mat] = build_S_matrices( ...
    uprior, uprior2, ...
    vprior, vprior2, ...
    wprior, wprior2, ...
    dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, ...
    dudx2, dudy2, dudz2, dvdx2, dvdy2, dvdz2, dwdx2, dwdy2, dwdz2, ...
    dudt, dvdt, dwdt, ...
    dudt2, dvdt2, dwdt2, ...
    is, ie, js, je, ks, ke)

    % Build s and s2 matrices similar to original code
    u_loc  = reshape(uprior(is:ie,js:je,ks:ke),  [], 1);
    v_loc  = reshape(vprior(is:ie,js:je,ks:ke),  [], 1);
    w_loc  = reshape(wprior(is:ie,js:je,ks:ke),  [], 1);
    u2_loc = reshape(uprior2(is:ie,js:je,ks:ke), [], 1);
    v2_loc = reshape(vprior2(is:ie,js:je,ks:ke), [], 1);
    w2_loc = reshape(wprior2(is:ie,js:je,ks:ke), [], 1);

    dudx_loc = reshape(dudx(is:ie,js:je,ks:ke), [], 1);
    dudy_loc = reshape(dudy(is:ie,js:je,ks:ke), [], 1);
    dudz_loc = reshape(dudz(is:ie,js:je,ks:ke), [], 1);
    dvdx_loc = reshape(dvdx(is:ie,js:je,ks:ke), [], 1);
    dvdy_loc = reshape(dvdy(is:ie,js:je,ks:ke), [], 1);
    dvdz_loc = reshape(dvdz(is:ie,js:je,ks:ke), [], 1);
    dwdx_loc = reshape(dwdx(is:ie,js:je,ks:ke), [], 1);
    dwdy_loc = reshape(dwdy(is:ie,js:je,ks:ke), [], 1);
    dwdz_loc = reshape(dwdz(is:ie,js:je,ks:ke), [], 1);

    dudx2_loc = reshape(dudx2(is:ie,js:je,ks:ke), [], 1);
    dudy2_loc = reshape(dudy2(is:ie,js:je,ks:ke), [], 1);
    dudz2_loc = reshape(dudz2(is:ie,js:je,ks:ke), [], 1);
    dvdx2_loc = reshape(dvdx2(is:ie,js:je,ks:ke), [], 1);
    dvdy2_loc = reshape(dvdy2(is:ie,js:je,ks:ke), [], 1);
    dvdz2_loc = reshape(dvdz2(is:ie,js:je,ks:ke), [], 1);
    dwdx2_loc = reshape(dwdx2(is:ie,js:je,ks:ke), [], 1);
    dwdy2_loc = reshape(dwdy2(is:ie,js:je,ks:ke), [], 1);
    dwdz2_loc = reshape(dwdz2(is:ie,js:je,ks:ke), [], 1);

    dudt_loc  = reshape(dudt(is:ie,js:je,ks:ke), [], 1);
    dvdt_loc  = reshape(dvdt(is:ie,js:je,ks:ke), [], 1);
    dwdt_loc  = reshape(dwdt(is:ie,js:je,ks:ke), [], 1);
    dudt2_loc = reshape(dudt2(is:ie,js:je,ks:ke), [], 1);
    dvdt2_loc = reshape(dvdt2(is:ie,js:je,ks:ke), [], 1);
    dwdt2_loc = reshape(dwdt2(is:ie,js:je,ks:ke), [], 1);

    n = numel(u_loc);
    s_mat  = zeros(n, 33);
    s2_mat = zeros(n, 15);

    s_mat(:,1:3)  = [u_loc, v_loc, w_loc];
    s2_mat(:,1:3) = [u2_loc, v2_loc, w2_loc];

    s_mat(:,4:12)  = [dudx_loc, dudy_loc, dudz_loc, ...
                      dvdx_loc, dvdy_loc, dvdz_loc, ...
                      dwdx_loc, dwdy_loc, dwdz_loc];

    s2_mat(:,4:12) = [dudx2_loc, dudy2_loc, dudz2_loc, ...
                      dvdx2_loc, dvdy2_loc, dvdz2_loc, ...
                      dwdx2_loc, dwdy2_loc, dwdz2_loc];

    s_mat(:,13:30) = 0;
    s_mat(:,31:33) = [dudt_loc, dvdt_loc, dwdt_loc];
    s2_mat(:,13:15)= [dudt2_loc, dvdt2_loc, dwdt2_loc];
end
