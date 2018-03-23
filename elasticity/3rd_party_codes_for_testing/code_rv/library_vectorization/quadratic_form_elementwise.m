function value_elementwise=quadratic_form_elementwise(u,elements,K_3Dmatrix)
    u_3Dvector=conv_ma2av(u(elements));
    value_elementwise=avtamav(u_3Dvector,K_3Dmatrix,u_3Dvector);
end
