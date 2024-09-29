function [output,omega]=real_valued_SBL(Y_bar,Fs,omega)
%%
[L,K]=size(Y_bar);
I=eye(L/2);
J=flip(I);
Q=[I,J;1j*J,-1j*I]/sqrt(2);
N=length(omega);
norm_y=norm(Y_bar,'fro')/sqrt(L*K);
Y_bar=Y_bar/norm_y;  
bar_A=exp(-1j*pi*(L-1:-2:1-L).'*omega/Fs)/sqrt(L);
bar_B=-1j*pi*(L-1:-2:1-L).'/Fs.*bar_A;
bar_A=real(Q*bar_A);
bar_B=real(Q*bar_B);
reslu=omega(2)-omega(1);
%% Initialization 
maxiter=200;
alpha=0.01;
gamma=ones(N,1);
a=1e-10;b=1e-10;
converged = false;
iter = 0;
etc=10;
It=150;
rho=0.9;
while ~converged && iter<maxiter   
    %% Equation (27)-(29)
    %% update q(X)     
    AA=bar_A.'*bar_A;
    Sigma=inv(alpha* AA  + diag(gamma));
    mu = Sigma * (alpha * (bar_A.' * Y_bar)   );  
    %% update q(alpha)  
    resid=Y_bar-bar_A*mu;
    term2=sum(diag( Sigma*AA));
    alpha_old=alpha;
    alpha=( L*K + 2*a )/( 2*b +  norm(resid, 'fro')^2+   K*real(term2) );
    alpha=rho*alpha_old+(1-rho)*alpha;
    %% update q(gamma) 
    Exx = sum( mu.*conj(mu),2) + K*real(diag(Sigma));
    c_k=K+2*a;
    d_k=real(2*b+Exx);
    gamma=c_k ./ d_k;
    %% Equation (34) 
    %% off-grid  
    if iter<It
        Pm=sum( mu.*conj(mu), 2);
        [~,sort_ind]=sort(Pm, 'descend');   
        idx=sort_ind(1:etc);
        BHB = bar_B(:,idx).' * bar_B(:,idx);
        P2= K * Sigma(idx,idx);  
        P = real( BHB .* ((mu(idx,:) * mu(idx,:).') +   P2   )  );
        v2= K * real(diag(bar_B(:,idx).' * bar_A * Sigma(:,idx)));  
        v = real(sum( mu(idx,:) .* (bar_B(:,idx).' * (Y_bar - bar_A * mu)),2)) -   v2;
        temp_grid=v./diag(P);
        temp_grid=temp_grid';
        theld=reslu/50*0.99^iter;
        temp_grid=sign(temp_grid)*theld;
        omega(idx)=omega(idx) + temp_grid;
        F_active=exp(-1j*pi*(L-1:-2:1-L).'*omega(idx)/Fs)/sqrt(L);
        bar_A(:,idx)=real(Q*F_active);
        bar_B(:,idx)=-1j*pi*(L-1:-2:1-L).'/Fs.*F_active;
        bar_B(:,idx)=real(Q*bar_B(:,idx));
    end
    %% Equation (36) 
    %%  reset the new grid
    if  iter==It   
        Pm=1./gamma;
        fn=omega(Pm==max(Pm(omega<6)));
        N=3;
        eta=1:1:N;
        fn_all= fn*eta;      
        omega=fn_all;
        A=exp(-1j*pi*(L-1:-2:1-L).'*omega/Fs)/sqrt(L);
        bar_A=real(Q*A);
        bar_B=real(Q*((-1j)*pi*(L-1:-2:1-L).'/Fs.*A));
        gamma=ones(length(omega),1);
        etc=length(omega);    
    end
    if iter>It
        rho=0.98;
        BHB = bar_B' * bar_B;
        P2= K * Sigma;  
        P = real( conj(BHB) .* ((mu * mu') +   P2   )  );
        v2= K * real(diag(bar_B' * bar_A * Sigma));  
        v = sum( real(conj(mu) .* (bar_B' * (Y_bar - bar_A * mu))),2) -   v2;
        temp_grid=(eta*v)/(eta*P*eta.');
        temp_grid=temp_grid';
        theld=reslu/50*0.99^(iter);
        temp_grid=sign(temp_grid)*theld;
        omega=omega + temp_grid*eta;
        F_active=exp(-1j*pi*(L-1:-2:1-L).'*omega/Fs)/sqrt(L);
        bar_A=real(Q*F_active);
        bar_B=real(Q*(((-1i)*pi*(L-1:-2:1-L).'/Fs).*F_active));
    end
    %% stopping criteria
    if iter >= maxiter
        converged = true;
    end
    iter = iter + 1;
end
output=sqrt(1./abs(gamma)*norm_y^2);  


