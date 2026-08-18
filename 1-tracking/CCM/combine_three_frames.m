function [Xp, Yp, Zp, Up, Vp, Wp, axp, ayp, azp, T, weight, ...
          UNCUp, UNCVp, UNCWp, UNCAxp, UNCAyp, UNCAzp] = ...
          combine_three_frames(last, now, next)

    % Concatenate three frames exactly as in original script.
    Xp  = [last.X;  now.X;  next.X];
    Yp  = [last.Y;  now.Y;  next.Y];
    Zp  = [last.Z;  now.Z;  next.Z];
    Up  = [last.U;  now.U;  next.U];
    Vp  = [last.V;  now.V;  next.V];
    Wp  = [last.W;  now.W;  next.W];
    axp = [last.ax; now.ax; next.ax];
    ayp = [last.ay; now.ay; next.ay];
    azp = [last.az; now.az; next.az];

    nLast = length(last.X);
    nNow  = length(now.X);
    nNext = length(next.X);

    T = [-ones(nLast,1); zeros(nNow,1); ones(nNext,1)];
    weight = ones(length(T),1);

    UNCUp  = [last.UNCU;  now.UNCU;  next.UNCU];
    UNCVp  = [last.UNCV;  now.UNCV;  next.UNCV];
    UNCWp  = [last.UNCW;  now.UNCW;  next.UNCW];
    UNCAxp = [last.UNCAx; now.UNCAx; next.UNCAx];
    UNCAyp = [last.UNCAy; now.UNCAy; next.UNCAy];
    UNCAzp = [last.UNCAz; now.UNCAz; next.UNCAz];
end
