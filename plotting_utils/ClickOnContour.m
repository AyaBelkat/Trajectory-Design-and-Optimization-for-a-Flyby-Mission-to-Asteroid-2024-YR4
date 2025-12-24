function ClickOnContour(~,~)

mu = 1.32712440018e11;  % Sun gravitational parameter


cp = get(gca,'CurrentPoint');   % coordinates of mouse click
xnum  = cp(1,1);           % clicked departure (datenum)
ynum  = cp(1,2);           % clicked arrival   (datenum)

% Convert to datetime (raw click, unsnapped)
dep_click_dt = datetime(xnum,'ConvertFrom','datenum');
arr_click_dt = datetime(ynum,'ConvertFrom','datenum');

fprintf('You clicked at (%s, %s)\n',dep_click_dt,arr_click_dt);


depEt = cspice_str2et(datestr(dep_click_dt));
arrEt = cspice_str2et(datestr(arr_click_dt));

D_earth = cspice_spkezr('Earth', depEt, 'J2000', 'NONE', 'Sun');
D_ast = cspice_spkezr('54509621',  arrEt, 'J2000', 'NONE', 'Sun');

TOF = arrEt - depEt;

r1 = D_earth(1:3);
r2 = D_ast(1:3);

vE = D_earth(4:6);
vA = D_ast(4:6);

[v1, v2] = lambert_solver(r1, r2, TOF, mu);


if size(v1, 1) == 1
    v1 = v1(:);
    v2 = v2(:);
    plot_lambert_trajectory(r1, r2, v1, 0);

elseif size(v1, 1) > 1
    v11 = v1(1, :);  v11 = v11(:);
    v21 = v2(1, :);  v21 = v21(:);

    v12 = v1(2, :);  v12 = v12(:);
    v22 = v2(2, :);  v22 = v22(:);

    v13 = v1(3, :);  v13 = v13(:);
    v23 = v2(3, :);  v23 = v23(:);

    
    C31 = norm(v11 - vE)^2;
    C32 = norm(v12 - vE)^2;
    C33 = norm(v13 - vE)^2;

    if C31 <= 10
        plot_lambert_trajectory(r1, r2, v11(:), 0);
        disp("C 0 rev is"); disp(C31);
        end
        if C32 <= 10
        plot_lambert_trajectory(r1, r2, v12(:), 1);
        disp("C 1 P rev is"); disp(C32);
        end
        if C33 <= 10
        plot_lambert_trajectory(r1, r2, v13(:), 1);
        disp("C 1 R rev is"); disp(C33);
    end
end





