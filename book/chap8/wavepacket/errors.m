
% Errors computed with different methods at time t=1 and t=2
%
%             dx          max-norm      1-norm

errmct1 = [ ...
            1.0000e-02   7.0796e-01   8.6617e-02 ; ...
	    5.0000e-03   2.5829e-01   2.2309e-02 ; ...
	    2.5000e-03   8.5244e-02   6.5538e-03 ; ...
	    1.2500e-03   2.8355e-02   1.9658e-03 ; ...
	    6.2500e-04   1.1062e-02   5.2749e-04 ; ...
	    3.1250e-04   4.6413e-03   1.3494e-04 ...
	  ];


errmct2 = [ ...
            1.0000e-02   8.5885e-01   1.0640e-01 ; ...
	    5.0000e-03   3.9795e-01   4.0521e-02 ; ...
	    2.5000e-03   1.2244e-01   9.6374e-03 ; ...
	    1.2500e-03   4.0759e-02   3.3628e-03 ; ...
	    6.2500e-04   1.6929e-02   9.7392e-04 ; ...
	    3.1250e-04   7.2366e-03   2.5607e-04 ...
	  ];


errlwt1 = [ ...
            1.0000e-02   1.0845e+00   1.4077e-01 ; ...
	    5.0000e-03   6.6908e-01   8.3756e-02 ; ...
	    2.5000e-03   2.0207e-01   2.4130e-02 ; ...
	    1.2500e-03   5.1820e-02   6.1444e-03 ; ...
	    6.2500e-04   1.3004e-02   1.5401e-03 ; ...
	    3.1250e-04   3.2533e-03   3.8517e-04  ...
	  ];


errlwt2 = [ ...
            1.0000e-02   9.9977e-01   1.2523e-01 ; ...
	    5.0000e-03   1.0485e+00   1.3710e-01 ; ...
	    2.5000e-03   3.8846e-01   4.7194e-02 ; ...
	    1.2500e-03   1.0318e-01   1.2264e-02 ; ...
	    6.2500e-04   2.5989e-02   3.0795e-03 ; ...
	    3.1250e-04   6.5059e-03   7.7032e-04  ...
	  ];

errsbt1 = [ ...
            1.0000e-02   6.2913e-01   7.5112e-02 ; ...
	    5.0000e-03   1.7050e-01   1.8688e-02 ; ...
	    2.5000e-03   6.6053e-02   5.0801e-03 ; ...
	    1.2500e-03   6.2229e-02   4.1482e-03 ; ...
	    6.2500e-04   3.6565e-02   1.8037e-03 ; ...
	    3.1250e-04   1.5949e-02   5.2307e-04  ...
	  ];


errsbt2 = [ ...
            1.0000e-02   8.0753e-01   9.8935e-02 ; ...
	    5.0000e-03   2.2674e-01   2.1296e-02 ; ...
	    2.5000e-03   8.2798e-02   7.1005e-03 ; ...
	    1.2500e-03   6.9484e-02   4.7224e-03 ; ...
	    6.2500e-04   5.4547e-02   3.0940e-03 ; ...
	    3.1250e-04   2.3720e-02   9.7326e-04 ...
	  ];

clf
%Haxes = axes('position',[.1 .1 .8 .8]);
set(gca,'fontsize',15)
loglog(errlwt1(:,1),errlwt1(:,2),'-','LineWidth',2.0)
hold on
loglog(errmct1(:,1),errmct1(:,2),'--','LineWidth',2.0)
%loglog(errsbt1(:,1),errsbt1(:,2),'-.')
title('max-norm errors at t = 1')
%legend('Lax-Wendroff','MC-limiter','Superbee')
legend('Lax-Wendroff','MC-limiter')
loglog(errlwt1(:,1),errlwt1(:,2),'.','MarkerSize',25)
loglog(errmct1(:,1),errmct1(:,2),'.','MarkerSize',25)

loglog([2e-4 2e-3],[2e-2 2e0])
text(2.3e-4,1e-1,'slope 2','fontsize',15)

hold off
query

clf
%Haxes = axes('position',[.1 .1 .8 .8]);
set(gca,'fontsize',15)
loglog(errlwt2(:,1),errlwt2(:,2),'-','LineWidth',2.0)
hold on
loglog(errmct2(:,1),errmct2(:,2),'--','LineWidth',2.0)
%loglog(errsbt2(:,1),errsbt2(:,2),'-.')
title('max-norm errors at t = 2')
%legend('Lax-Wendroff','MC-limiter','Superbee')
legend('Lax-Wendroff','MC-limiter')
loglog(errlwt2(:,1),errlwt2(:,2),'.','MarkerSize',25)
loglog(errmct2(:,1),errmct2(:,2),'.','MarkerSize',25)

loglog([2e-4 2e-3],[2e-2 2e0])
text(2.3e-4,1e-1,'slope 2','fontsize',15)

hold off

ne = size(errmct2,1);
A = [1  log(errmct2(ne-1,1));  1 log(errmct2(ne,1))];
b = [log(errmct2(ne-1,2));  log(errmct2(ne,2))];
y = A\b;
C = exp(y(1));
order = y(2);
disp(['max-norm error for MC = ' num2str(C) ' * h^(' num2str(order) ')'])

ne = size(errlwt2,1);
A = [1  log(errlwt2(ne-1,1));  1 log(errlwt2(ne,1))];
b = [log(errlwt2(ne-1,2));  log(errlwt2(ne,2))];
y = A\b;
C = exp(y(1));
order = y(2);
disp(['max-norm error for LW = ' num2str(C) ' * h^(' num2str(order) ')'])

query

clf
%Haxes = axes('position',[.1 .1 .8 .8]);
set(gca,'fontsize',15)
loglog(errlwt1(:,1),errlwt1(:,3),'-','LineWidth',2.0)
hold on
loglog(errmct1(:,1),errmct1(:,3),'--','LineWidth',2.0)
%loglog(errsbt1(:,1),errsbt1(:,3),'-.')
title('1-norm errors at t = 1')
%legend('Lax-Wendroff','MC-limiter','Superbee')
legend('Lax-Wendroff','MC-limiter')
loglog(errlwt1(:,1),errlwt1(:,3),'.','MarkerSize',25)
loglog(errmct1(:,1),errmct1(:,3),'.','MarkerSize',25)

loglog([2e-4 2e-3],[2e-2 2e0])
text(2.3e-4,1e-1,'slope 2','fontsize',15)

hold off
query

clf
%Haxes = axes('position',[.1 .1 .8 .8]);
set(gca,'fontsize',15)
loglog(errlwt2(:,1),errlwt2(:,3),'-','LineWidth',2.0)
hold on
loglog(errmct2(:,1),errmct2(:,3),'--','LineWidth',2.0)
%loglog(errsbt2(:,1),errsbt2(:,3),'-.')
title('1-norm errors at t = 2')
%legend('Lax-Wendroff','MC-limiter','Superbee')
legend('Lax-Wendroff','MC-limiter')
loglog(errlwt2(:,1),errlwt2(:,3),'.','MarkerSize',25)
loglog(errmct2(:,1),errmct2(:,3),'.','MarkerSize',25)

loglog([2e-4 2e-3],[2e-2 2e0])
text(2.3e-4,1e-1,'slope 2','fontsize',15)

hold off

ne = size(errmct2,1);
A = [1  log(errmct2(ne-1,1));  1 log(errmct2(ne,1))];
b = [log(errmct2(ne-1,3));  log(errmct2(ne,3))];
y = A\b;
C = exp(y(1));
order = y(2);
disp(['1-norm error for MC = ' num2str(C) ' * h^(' num2str(order) ')'])

ne = size(errlwt2,1);
A = [1  log(errlwt2(ne-1,1));  1 log(errlwt2(ne,1))];
b = [log(errlwt2(ne-1,3));  log(errlwt2(ne,3))];
y = A\b;
C = exp(y(1));
order = y(2);
disp(['1-norm error for LW = ' num2str(C) ' * h^(' num2str(order) ')'])