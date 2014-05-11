%Belief Propagation in Log Space

%Inputs:
%E: an M by 2 list of edges (M is the number of edges). Edges should be
%redundant, if (s,t) is in the list, so is (t,s).
%U: an N by q array providing univariate functions
%B: an M by q by q array providing bivariate functions
%mode: mode of operation. Possible values are 1 (fing marginal
%probabilities), 2 (find-max marginal probabilities), 3 (find normalizing
%constant).

%Outputs:
%Prbs: the N by q list of associated probabilities

function Prbs = beliefPropLog(E,U,B,mode)

N = size(U,1); %number of nodes
M = size(E,1); %number of edges
q = size(U,2); %number of states

%depending on mode, setup up appropriate function handles
if mode == 1,
    sumOp = @logsumexp;
elseif mode == 2,
    sumOp = @max;
elseif mode == 3,
    sumOp = @logsumexp;
end;

%initialize messages
Msgs = ones(M,q)/q;

%loop over edges
epsilon = 1e-10; %determines when to stop
changes = [];
for t = 1:100,
    old = Msgs;
    for i = 1:M,
        
        %find incoming edges
        from = E(i,1);
        to = E(i,2);
        in = find(E(:,1)~=to & E(:,2)==from);
        %pass message
        for j = 1:q,
            msgs = [0 0];
            for k = 1:q,
                msgs(k) = log(B(i,k,j))+log(U(from,k))+sum(Msgs(in,k));
            end;
            Msgs(i,j) = sumOp(msgs);
        end;        
    end;
    change = sum(sum(abs(Msgs-old)));
    changes = [changes change];
    if change < epsilon, break; end;
end;
%check convergence
%plot(changes)
%figure;
%compute probabilities
Prbs = zeros(N,q);
for i = 1:N,
    in = find(E(:,2)==i);
    for j = 1:q,
        Prbs(i,j) = exp(log(U(i,j))+sum(Msgs(in,j)));
    end;
    %normalize
    if mode ~= 3,
        Prbs(i,:) = Prbs(i,:)/sum(Prbs(i,:));
    end;
end;

%return normalizing constant Z if that's what was asked for
if mode == 3,
    Prbs = sum(Prbs,2);
    Prbs = Prbs(1);
end;

%numerically stable log of sum of exponentials
function y = logsumexp(x)

m = max(x);
x = x-m;
y = m + log(sum(exp(x)));



