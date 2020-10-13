function [X,D,e] = sparse_solver_irls(p, A, Y, beta, terminationValue, maxIterations)
%SPARSE_SOLVER_IRLS Steered-response power map using a MVDR beamformer
%   
%   This routine performs sparse recovery using the iterative reweighted
%   least squares (IRLS) algorithm, with an L_{p,2} norm, assuming
%   spatiotemporal signals Y, of M channels/sensors/microphones and T
%   temporal frames (snapshots), and a linear mixing model of the form 
%   Y = A*X, where A is a MxK overcomplete basis dictionary, with K>>M. 
%   The original algorithm is described in:
%   
%   Chartrand, R., & Yin, W. (2008). Iteratively reweighted algorithms for 
%   compressive sensing. In 2008 IEEE International Conference on Acoustics, 
%   Speech and Signal Processing (ICASSP), pp. 3869-3872.
%
%   and further inspired by the microphone array processing work of
%
%   Wabnitz, A., Epain, N., McEwan, A., & Jin, C. (2011). Upscaling ambisonic 
%   sound scenes using compressed sensing techniques. In 2011 IEEE Workshop 
%   on Applications of Signal Processing to Audio and Acoustics (WASPAA).
%
%   Inputs:
%       p:  sparsity exponent p<=1, of the L_{p,2} norm
%       A:  overcomplete dictionary (in CS jargon) of steering vectors for
%           a dense grid of K directions, with A of size MxK.
%       Y:  measurement vectors (in CS jargon), multichannel microphone signals 
%           in audio/acoustics, of size MxT
%       beta:   regularization parameter
%       terminationValue:   error threshold to stop iterations in the
%                           solver
%       maxIterations:  maximum number of iterations before stopping the
%                       solver, if terminationValue is not reached
%
%   Outputs:
%       X:  KxT sparse signal solution
%       D:  demixing matrix of size KxM
%       e:  Kx1 sparse signal powers
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SPARSE_SOLVER_IRLS.M - 13/5/2019
% Archontis Politis, archontis.politis@tuni.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<6
    maxIterations = 0;
end

nChan = size(A,1); % number of channels/sensors
nDict = size(A,2); % size of steering vector dictionary
nSnap = size(Y,2); % size of temporal snapshots (e.g. STFT frames, or time-domain samples)
% initial values for amplitudes and weights
W = eye(nDict);
minEpsilon = terminationValue;
epsilon = 1;
n = 0;

if nSnap>nChan
    % dimensionality reduction using SVD
    [U, L, ~] = svd(Y);
    L_trunc = L(:,1:nChan);
    Y2 = U*L_trunc;
else
    Y2 = Y;
end
X_prev = A\Y2; % initial LS solution
while epsilon>minEpsilon
    % regularization
    gamma = beta/(1-beta) * trace(A*W*A')/nChan;
    % gamma = beta^0.15/(1-beta^0.15) * trace(A*W*A')/nChan; % Epain & Jin
    % reconstruction matrix
    D = (W*A')/(A*W*A' + gamma*eye(nChan));
    % amplitude estimates
    X_current = D * Y2;    
    % norm of difference between previous and current solution
    er_k = sqrt(sum((X_current - X_prev).^2,2));
    er = sqrt(sum(er_k.^2));
    % sparse solution signal energies
    e = sum(X_current.^2,2);
    e_max = max(e);
    epsilon = e_max/nDict; % Epain & Jin
%     if er < epsilon^(1/2)/100 % Chartrand & Yin
%         epsilon = epsilon /10;
%     end
    w = (e + epsilon).^(1-p/2);
    W = diag(w);
    
    n = n+1;
    disp(['Iteration ' num2str(n) ': ' num2str(epsilon)])
    if maxIterations
        if n>maxIterations, break; end
    end
end
if nSnap>nChan
    D = X_current*inv(U*L_trunc);
    X = D*Y;
else
    X = X_current;
end
