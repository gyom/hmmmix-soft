function [gamma,alpha,beta,loglik] = fwd_back(init_state_distrib, transmat, obslik)

% forwards propagation, backwards sampling
%
% input
% init_state_distrib(q)
% transmat(q, q')
% obslik(q, t)
% nsamples
%
% output
% samples(t, s)  = sample s of discrete state at time  t


[Q T] = size(obslik);
scale = ones(1,T);
loglik = 0;
alpha = zeros(Q,T, 'single');
beta = zeros(Q,T,'single');
gamma = zeros(Q,T,'single');
trans = transmat;

t = 1;
alpha(:,1) = init_state_distrib(:) .* obslik(:,t);
[alpha(:,t), scale(t)] = normalize(alpha(:,t));
for t=2:T
  m = trans' * alpha(:,t-1);
  alpha(:,t) = m(:) .* obslik(:,t);
  [alpha(:,t), scale(t)] = normalize(alpha(:,t));
  %assert(approxeq(sum(alpha(:,t)),1))
end
loglik = sum(log(scale));


beta = zeros(Q,T);
t=T;
beta(:,T) = ones(Q,1);
gamma(:,T) = alpha(:,T);

for t=T-1:-1:1
 b = beta(:,t+1) .* obslik(:,t+1);
 beta(:,t) = normalize(transmat * b);
 gamma(:,t) = normalize(alpha(:,t) .* beta(:,t));
end
