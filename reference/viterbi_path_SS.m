function [path loglik seg] = viterbi_path_SS(prior, transmat, obslik)
% VITERBI Find the most-probable (Viterbi) path through the HMM state trellis.
% path = viterbi(prior, transmat, obslik)
%
% Inputs:
% prior(i) = Pr(Q(1) = i)
% transmat(i,j) = Pr(Q(t+1)=j | Q(t)=i)
% obslik(i,t) = Pr(y(t) | Q(t)=i)
%
% Outputs:
% path(t) = q(t), where q1 ... qT is the argmax of the above expression.


% delta(j,t) = prob. of the best sequence of length t-1 and then going to state j, and O(1:t)
% psi(j,t) = the best predecessor state, given that we ended up in state j at t

scaled = 0;

T = size(obslik, 2);
prior = prior(:);
Q = length(prior);

delta = zeros(Q,T);
psi = zeros(Q,T);
path = zeros(1,T);
scale = ones(1,T);


t=1;
delta(:,t) = log(prior) + log(obslik(:,t));
if scaled
  [delta(:,t), n] = normalise(delta(:,t));
  scale(t) = 1/n;
end
psi(:,t) = 0; % arbitrary value, since there is no predecessor to t=1
for t=2:T
  for j=1:Q
    [delta(j,t), psi(j,t)] = max(delta(:,t-1) + log(transmat(:,j)));
    delta(j,t) = delta(j,t) + log(obslik(j,t));
  end
  if scaled
    [delta(:,t), n] = normalise(delta(:,t));
    scale(t) = 1/n;
  end
end
[p, path(T)] = max(delta(:,T));
seg = [1,T,path(T),0];
for t=T-1:-1:1
  path(t) = psi(path(t+1),t+1);
  if path(t) ~= path(t+1)
      seg(1,1)=t+1;
      seg(1,4) = 0; %computeSegmentBayesFactor(path(t+1),transmat,obslik(:,seg(end,1):seg(end,2)));
      seg = [[1,t,path(t),0];seg];
  end
end

% If scaled==0, p = prob_path(best_path)
% If scaled==1, p = Pr(replace sum with max and proceed as in the scaled forwards algo)
% Both are different from p(data) as computed using the sum-product (forwards) algorithm


if scaled
  loglik = -sum(log(scale));
  %loglik = prob_path(prior, transmat, obslik, path);
else
  loglik = p;
end
