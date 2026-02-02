function results=run_solver(opts)
    %RUN_SOLVER  Core solver driver (Algorithm 3) wrapped as a function.
    %
    % Input:
    %   opts : struct of parameters. Expected fields include:
    %          N, r_l, m, mu, tau, maxiter, maxiter_opt, sigma
    %
    % Output:
    %   results : full path of the saved .mat result file,
    %             produced by result_filename(opts).

    % -------------------------
    % 1) Unpack options 
    % -------------------------
    N=opts.N;
    r_l=opts.r_l;
    m=opts.m;
    mu=opts.mu;
    tau=opts.tau;
    maxiter=opts.maxiter;
    maxiter_opt=opts.maxiter_opt;
    sigma=opts.sigma;
    
    % -------------------------
    % 2) Load pre-generated data and initialization
    % -------------------------
    load(data_filename(opts));
    load(s0_filename(opts));
    
    % -------------------------
    % 3) Build Manopt manifold and problem container
    % -------------------------
    % Variable yt consists of:
    %   yt.y : complex Euclidean variable of size (3N) x (r_l*m)
    %   yt.t : complex Euclidean variable of size (3N+1) x 1
    elems.y=euclideancomplexfactory(3*N, r_l*m);
    elems.t=euclideancomplexfactory(3*N+1,1);
    manifold=productmanifold(elems);
    problem.M=manifold;
    
    % Manopt solver options for inner loop (rlbfgs)
    options.maxiter=maxiter_opt;
    options.memory=5;
    
    % -------------------------
    % 4) Outer-loop initialization: compute C, X, update M, and compute initial residuals
    % -------------------------
    C=J-H(yt.y,yt.t,N,m,r_l);
    
    X=R*C;
    
    % Update M using current sigma 
    M=M-sigma*X;
    
    M_mat=reshape(M,3*N+1,3*N+1);
    
    % Compute largest/smallest eigenvalues needed by eta_P
    [Vmax,lambda_max] = eigs(M_mat, 1, 'LA', eigopts);
    [Vmin,lambda_min] = eigs(M_mat, 1, 'SA', eigopts);
    
    iter=1;
    
    % eta_P: primal feasibility
    eta_P(iter)=max(0,-real(lambda_min))/(1+max(0,real(lambda_max)));
    
    % eta_D: dual feasibility
    eta_D(iter)=norm(X)*cD;
    
    % primal_obj: primal objective
    primal_obj=real(J'*M);
    
    % dual_obj: dual objective
    dual_obj=sum(C(idx_u));
    
    % eta_g: duality gap 
    eta_g(iter)=abs(primal_obj-dual_obj)/(1+abs(primal_obj)+abs(dual_obj));
    
    energy(iter)=primal_obj;
    
    tic
    
    % -------------------------
    % 5) Outer loop: ALM updates with Manopt inner solve
    % -------------------------
    % Stop when all three measures <= 1e-3 or when iter exceeds maxiter.
    while max([eta_P(iter),eta_D(iter),eta_g(iter)])>1e-3 && iter<=maxiter
        iter=iter+1;
        
        % Define Manopt cost and Euclidean gradient at current sigma/M
        % Loss and dLoss both accept a 'store' argument to enable caching.
        problem.cost=@(yt,store) Loss(yt.y,yt.t, N, m, r_l, sigma, J, M, R, idx_u, store);
        problem.egrad = @(yt,store) dLoss(yt.y,yt.t, N, m, r_l, sigma, J, M, R, u, V, W, store);
    
        % Inner solve: update yt to approximately minimize Loss
        [yt,~]=rlbfgs(problem,yt,options);
     
        C=J-H(yt.y,yt.t,N,m,r_l);
        X=R*C;
        
        % Update M
        M=M-sigma*X;
    
        M_mat=reshape(M,3*N+1,3*N+1);
       
        % Warm-start eigs with previous eigenvectors for speed
        eigopts.v0 = Vmax;
        [Vmax,lambda_max] = eigs(M_mat, 1, 'LA', eigopts);
        eigopts.v0 = Vmin;
        [Vmin,lambda_min] = eigs(M_mat, 1, 'SA', eigopts);
    
        % Update feasibility measures and objective values
        eta_P(iter)=max(0,-real(lambda_min))/(1+max(0,real(lambda_max)));
        eta_D(iter)=norm(X)*cD;
        primal_obj=real(J'*M);
        dual_obj=sum(C(idx_u));
        eta_g(iter)=abs(primal_obj-dual_obj)/(1+abs(primal_obj)+abs(dual_obj));
        energy(iter)=primal_obj;
    
        % -------------------------
        % 6) Adaptive update of sigma based on ratio rho = eta_P / eta_D
        % (residual balance strategy)
        % -------------------------
        rho=eta_P(iter)/eta_D(iter);
        if rho > mu 
            sigma=sigma/tau;
        elseif rho < 1/mu 
            sigma=sigma*tau;
        end
        
        % Print progress
        time=toc;
        fprintf('Iteration: %f\n',iter);
        fprintf('Primal feasibility measure: %f\n',log10(eta_P(iter)));
        fprintf('Dual feasibility measure: %f\n',log10(eta_D(iter)));
        fprintf('Duality Gap: %f\n',log10(eta_g(iter)));
        fprintf('Time: %f\n',time);
    end
    
    % Total time
    time=toc;
    
    % -------------------------
    % 7) Post-run summary statistics
    % -------------------------
    % relative_err 
    relative_err=abs(true_E0+energy(iter))/true_E0;
    
    % Print final summary
    fprintf('Time: %f\n',time);
    fprintf('Iteration: %f\n',iter);
    fprintf('Estimated Ground Energy : %f\n',primal_obj);
    fprintf('Relative Error: %f\n',relative_err);
    fprintf('Primal Feasibility Measure: %f\n',log10(eta_P(iter)));
    fprintf('Dual Feasibility Measure: %f\n',log10(eta_D(iter)));
    fprintf('Duality Gap: %f\n',log10(eta_g(iter)));
    
    % energy_per_site: magnitude of per-iteration energy change normalized by N
    energy_per_site = zeros(maxiter,1);
    energy_per_site(2:iter) = abs(diff(energy(1:iter))) / N;
    
    % -------------------------
    % 8) Plot diagnostics (log10 scale)
    % -------------------------
    close all;
    figure; hold on;
    plot(log10(eta_P))
    plot(log10(eta_D))
    plot(log10(eta_g))
    plot(log10(energy_per_site))
    hold off;
    
    % -------------------------
    % 9) Save results to disk
    % -------------------------
    results=result_filename(opts);
    save(results, "eta_P", "eta_D", "eta_g", "energy", "energy_per_site", ...
        "r_l","tau","maxiter_opt","mu","iter","relative_err","time");
end
