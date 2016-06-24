function exact_soln = exact3(x,time, K,m,x_ori)
% wave_speed = feval(speed{n_speed},x,time);
hackFactor=2;
% if m==0.5
%     char_curve = x*exp(K*time);
% else
    char_curve =(x.^(1-(hackFactor*m))+(1-(hackFactor*m))*K*time).^(1/(1-(hackFactor*m)));
% end
% error(num2str(time))
exact_soln = feval(@shape3,char_curve,x_ori);
return;