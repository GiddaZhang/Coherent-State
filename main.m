beta = 0.2;
alpha = 2;
w = 1;
N = 30;
S1 = 100;
S2 = 50;
XLim = 25;
psi = zeros(N, S1);
coff = zeros(N, 1);
x = linspace(-XLim, XLim, S1)';
for n = 1:N
    [psi(n,:),x] = get_state(n-1, beta, -XLim, XLim, S1);
    coff(n) = alpha^(n-1)/sqrt(factorial(n-1))*exp(-1/2*abs(alpha)^2);
end

t = linspace(0, 2*pi/w, S2);
e = zeros(N, S2);
for n = 1:N
    e(n,:) = coff(n)*exp(-1j*(n-1+1/2)*w*t);
end
base_line = zeros(1, S1);
for n = 1:N
    base_line = base_line + abs(coff(n))^2*psi(n,:).^2;
end
figure;
for t_ = 1:length(t)
    psi_abs = zeros(1, S1);
    for n = 1:N
        psi_abs = psi_abs + psi(n,:)*e(n, t_);
    end
    psi_abs = abs(psi_abs).^2;

    plot(x, psi_abs, x, base_line, 'LineWidth', 1.5);
    ylim([0 0.15]);
    legend('$|\psi|^2$','$\sum_n |c_n|^2\psi_n^2$', 'Interpreter', 'latex');
    pause(0.1);
end

function [psi_n, x] = get_state ...
    (n, alpha, x_l, x_r, sample)
    x = linspace(x_l, x_r, sample)';
    N_n = sqrt(alpha/sqrt(pi)/2^n/factorial(n));
    H_n = zeros(sample, 1);
    for idx = 1:sample
        H_n(idx) = hermite(n, alpha*x(idx));
    end
    psi_n = H_n*N_n.*exp(-1/2*alpha^2*x.^2);
end

function H = hermite(n, x)
    syms xi
    d = diff(exp(-xi^2), n);
    d_ = subs(d, 'xi', x);
    H = (-1)^n*exp(x^2)*d_;
end