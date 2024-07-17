function [x_best,f_best,counter,f_conv_history,previous_iters,previous_f]=particleSwarmWithRestarts(func,X_0_best,X_lb,X_ub,N_particles,max_iter,func_tol,x_tol,weight_mat,phi_g,phi_p,existing_file,save_file)
%% Particle Swarm Optimization function with built in saving and loading from previous restart
% Written by Thor E. Andreassen, PhD
%% Last Edited 7/16/24
disp('test')
if nargin<=6
    func_tol=.001;
    x_tol=.001;
end
if nargin <=8
    weight_mat=.1*ones(length(X_0_best),1);
    phi_g=.01;
    phi_p=.01;

end

if nargin <12
    existing_file='';
    save_file='';
end
if ~isempty(existing_file)
    try
        func_tol2=func_tol;
        x_tol2=x_tol;
        load(existing_file);
        func_tol=func_tol2;
        x_tol=x_tol2;
        clear func_tol2;
        clear x_tol2;
        loaded_data=1;
    catch
        loaded_data=0;
        try
            clear func_tol2;
            clear x_tol2;
        catch
            %
        end
    end
else
    loaded_data=0;
end

% create initial random starting locations
if loaded_data==0
    for count_particle=1:N_particles
        for count_designvar=1:length(X_ub)
            current_rand_pos=rand.*(X_ub(count_designvar)-X_lb(count_designvar))+X_lb(count_designvar);
            X_particle_curr(count_particle,count_designvar)=current_rand_pos;

            current_rand_vel=2*rand*(norm(X_ub(count_designvar)-X_lb(count_designvar)))-...
                norm(X_ub(count_designvar)-X_lb(count_designvar));
            V_particle_curr(count_particle,count_designvar)=current_rand_vel;
        end
    end

    % lh_val=lhsdesign(N_particles,length(X_ub));
    %
    % for count_particle=1:N_particles
    %     for count_designvar=1:length(X_ub)
    %         current_rand_pos=lh_val(count_particle,count_designvar).*(X_ub(count_designvar)-X_lb(count_designvar))+X_lb(count_designvar);
    %         X_particle_curr(count_particle,count_designvar)=current_rand_pos;
    %     end
    % end



    X_particle_best=X_particle_curr;
    X_particle_best_f=ones(N_particles,1)*Inf;
    X_global_best=X_0_best;
    X_global_best_f=func(X_0_best');
    X_global_best_f_history=[];
    previous_iters=[];
    previous_f=[];

    counter=1;
    func_conv=Inf;
    x_conv=Inf;
end


alpha=0.25;

if nargin>6
    while counter<=max_iter && func_conv>func_tol && x_conv>x_tol
        for count_particle=1:N_particles
            if counter~=1
                for count_designvar=1:length(X_ub)
                    rp=rand;
                    rg=rand;
                    V_particle_curr(count_particle,count_designvar)=...
                        weight_mat(count_designvar)*V_particle_curr(count_particle,count_designvar)+...
                        phi_p*rp*(X_particle_best(count_particle,count_designvar)-X_particle_curr(count_particle,count_designvar))+...
                        phi_g*rg*(X_global_best(count_designvar)-X_particle_curr(count_particle,count_designvar));
                    X_particle_curr(count_particle,count_designvar)=coerce(X_particle_curr(count_particle,count_designvar)+V_particle_curr(count_particle,count_designvar)*alpha,X_lb(count_designvar),X_ub(count_designvar));


                end
            end

            current_f=func(X_particle_curr(count_particle,:));
            previous_iters=[previous_iters,X_particle_curr(count_particle,:)'];
            previous_f=[previous_f,current_f];
            if ~isempty(save_file)
                save(save_file);
            end

            % current particle better
            if current_f < X_particle_best_f(count_particle)
                X_particle_best_f(count_particle)=current_f;
                X_particle_best(count_particle,:)=X_particle_curr(count_particle,:);
            end

            % current global better
            if current_f < X_global_best_f
                func_conv=abs(current_f-X_global_best_f)
                X_global_best_f=current_f;
                x_conv=norm(reshape(X_particle_curr(count_particle,:),size(X_global_best,1),size(X_global_best,2))-X_global_best);
                X_global_best(:)=reshape(X_particle_curr(count_particle,:),size(X_global_best,1),size(X_global_best,2));
            end
        end
        counter=counter+1;
        X_global_best_f_history=[X_global_best_f_history,X_global_best_f];

    end
    x_best=X_global_best;
    f_best=X_global_best_f;
    f_conv_history=X_global_best_f_history;
else
    while counter<=max_iter
        for count_particle=1:N_particles
            for count_designvar=1:length(X_ub)
                rp=rand;
                rg=rand;
                V_particle_curr(count_particle,count_designvar)=...
                    weight_mat(count_designvar)*V_particle_curr(count_particle,count_designvar)+...
                    phi_p*rp*(X_particle_best(count_particle,count_designvar)-X_particle_curr(count_particle,count_designvar))+...
                    phi_g*rg*(X_global_best(count_designvar)-X_particle_curr(count_particle,count_designvar));
                X_particle_curr(count_particle,count_designvar)=coerce(X_particle_curr(count_particle,count_designvar)+V_particle_curr(count_particle,count_designvar)*alpha,X_lb(count_designvar),X_ub(count_designvar));

            end

            current_f=func(X_particle_curr(count_particle,:));
            previous_iters=[previous_iters,X_particle_curr(count_particle,:)'];
            previous_f=[previous_f,current_f];
            if ~isempty(save_file)
                save(save_file);
            end

            % current particle better
            if current_f < X_particle_best_f(count_particle)
                X_particle_best_f(count_particle)=current_f;
                X_particle_best(count_particle,:)=X_particle_curr(count_particle,:);
            end

            % current global better
            if current_f < X_global_best_f
                func_conv=abs(current_f-X_global_best_f)
                X_global_best_f=current_f;
                x_conv=norm(reshape(X_particle_curr(count_particle,:),size(X_global_best,1),size(X_global_best,2))-X_global_best);
                X_global_best(:)=reshape(X_particle_curr(count_particle,:),size(X_global_best,1),size(X_global_best,2));
            end
        end
        counter=counter+1;
        X_global_best_f_history=[X_global_best_f_history,X_global_best_f];

    end
    x_best=X_global_best;
    f_best=X_global_best_f;
    f_conv_history=X_global_best_f_history;
end
end

function y = coerce(x,bl,bu)
% return bounded value clipped between bl and bu
y=min(max(x,bl),bu);
end