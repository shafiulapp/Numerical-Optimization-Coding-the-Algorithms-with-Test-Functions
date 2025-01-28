classdef algorithms < handle
  
    properties     
        plot_sum = plotsum;           
        num_step = 0;        
        num_trs = 0;              
        f ;        
        plots = true;        
        opt_x; 
        max_iter = 10000; 
    end
    
    methods
        function run(obj, x_0, f_kind, method, plots)           
            obj.plots = plots;           
            obj.f = f_obj(f_kind);    
            obj.plot_sum = plotsum();
            
            switch (method)
                case "SDLS"       
                    tic
                    SD_LS(obj, x_0);
                case "CG"          
                    tic
                    CG(obj, x_0, "FR"); %PR %FR
                case "NewtonLS"   
                    tic
                    Newton_LS(obj, x_0);
                case "TRS"    
                    tic
                    Trust_Region(obj, x_0,"iterative"); %dogleg %iterative
                case "QN"          
                    tic
                    Quasinewton(obj, x_0, "DFP"); %BFGS %DFP
                case "IQN"         
                    tic
                    inexactnewton(obj, x_0);
                otherwise
                    disp("Incorrect method name");
            end
            toc

            disp("# function evaluations: " + obj.f.fval_count);
           disp("# function gradient evaluations: " + obj.f.gval_count);
            disp("# function hessian evaluations: " + obj.f.hval_count);
            disp("# iterations: " + obj.num_step);
           disp("# trust region subproblems: " + obj.num_trs);
           disp("Optimal x: [" + num2str(obj.opt_x(1)) + ", " + num2str(obj.opt_x(2)) + "]");
            
            if (obj.plots)
                obj.plot_sum();
            end
        end
        
     
        
       

       
        
        
        function SD_LS(obj, x_0)           
            x = x_0;
            fg = obj.f.grad(x);
            a = 1;
            p = [2,2].';  % initial direction
            
            while (criteria_stop(obj, fg, a*p))
               if obj.num_step > obj.max_iter
                 break;
               end
                fg = obj.f.grad(x);
                p = -fg;  %desecent direction
                a = linesearch(obj, p, x,"wolfe",true);
               
                if (obj.plots)
                    obj.plot_sum.fvals(end+1) = obj.f.val(x);
                    obj.plot_sum.gnorm(end+1) = norm(fg);
                    obj.plot_sum.step_lengths(end+1) = norm(a*p);
                end
       
                x = x + a*p;
                obj.opt_x = x;      
                obj.num_step = obj.num_step + 1;
            end
        end
        


%algorithm 5.4 pg 121
        function CG(obj, x_0, submethod)
            
            x = x_0;
            fg = obj.f.grad(x);
            p = -fg;
            a = 1;
            
            while (criteria_stop(obj, fg, a*p))
                if obj.num_step > obj.max_iter
                 break;
                end
                a = linesearch(obj, p, x,"wolfe",true);% ref: pg 122 to use wolfe
                x = x + a*p;
                f2_g = obj.f.grad(x);
                beta = 0;
                if (submethod == "FR")
                    beta = (f2_g.' * f2_g) / (fg.' * fg); %algorithm 5.4 page 121 eq 5.41a
                elseif (submethod == "PR")
                    beta = (f2_g.' * (f2_g - fg)) / (fg.' * fg); %eq 5.44 pg 122
                end
              
                if (obj.plots)
                    obj.plot_sum.fvals(end+1) = obj.f.val(x);
                    obj.plot_sum.gnorm(end+1) = norm(fg);
                    obj.plot_sum.step_lengths(end+1) = norm(a*p);
                end                                
                p = beta*p - f2_g;
                fg = f2_g;
                obj.opt_x = x;  % Store the current solution
                obj.num_step = obj.num_step + 1;
            end
        end






        function Newton_LS(obj, x_0)        
            x = x_0;
            fg = obj.f.grad(x);
            a = 1;
           p = [2,2].';  
            
            while (criteria_stop(obj, fg, a*p))
                if obj.num_step > obj.max_iter
                 break;
                end
                fg = obj.f.grad(x);
                B = obj.f.hess(x);  % exact Hessian

                p = -B \ fg;
                a = linesearch(obj, p, x,"wolfe",true);

                if (obj.plots)
                    obj.plot_sum.fvals(end+1) = obj.f.val(x);
                    obj.plot_sum.gnorm(end+1) = norm(fg);
                    obj.plot_sum.step_lengths(end+1) = norm(a*p);
                end

                x = x + a*p;
                obj.opt_x = x;        
                obj.num_step = obj.num_step + 1;
            end
        end
      


%alg 4.1
               function Trust_Region(obj,x_0,submethod)
            x=x_0;
            p=[2,2].';% need to check
            eta=0.1;
            del_hat=2;
            delta=del_hat*0.5;        
            fg=obj.f.grad(x);
            a=2;
             while(criteria_stop(obj,fg, a*p))
                 if obj.num_step > obj.max_iter
                 break;
                 end
                obj.num_step=obj.num_step+1;
                g=obj.f.grad(x);
                H=obj.f.hess(x);
                m_y=@(y) obj.f.val(x)+(g.'*y)+(1/2)*y.'*H*y;



                if submethod=="dogleg"
                    pu=-((g.'*g)/(g.'*H*g))*g; %eq 4.15
                    pb=-H\g;
                    p=piecewise_dogleg(obj,pb,pu,delta);

                elseif submethod=="iterative"
                    lamda_0=-min(eig(H))+1e-6; %initial guess
                    lamda=iteratives(obj,lamda_0,delta,H,g);
                    p=(H+(lamda*eye(2)))\(-g); %eq4.37
                end





                % 
                %  try
                % 
                %     pu = -((g.' * g) / (g.' * H * g)) * g;  %eq 4.15
                %     pb = -H \ g; % Newton direction
                % 
                % 
                % 
                %     if all(eig(H)> 0)
                %         p = piecewise_dogleg(obj, pb, pu, delta);
                %     else
                %         error('Hessian is not positive definite, switching to iterative method.');
                %     end
                % catch
                %     % Fallback to iterative method 
                %     lamda_0 = -min(eig(H)) + 0.1; % Initial guess 
                %     lamda = iteratives(obj, lamda_0, delta, H, g);
                %     p = (H + (lamda * eye(2))) \ (-g); %before eq 4.37
                %  end



                obj.num_trs=obj.num_trs+1;
                if p==0
                    rho_k=0;
                else 
                    f_x=obj.f.val(x);
                    f_p=obj.f.val(x+p);
                    m_0=m_y([0,0].');
                    m_p=m_y(p);
                    rho_k=(f_x-f_p)/(m_0-m_p); %ratio actual/pred
                end

                if rho_k<0.25
                    delta=delta*0.25;
                else
                    if rho_k>0.75 && norm(p)>delta 
                        delta=min(2*delta,del_hat);
                    end
                end
                if rho_k>eta
                    x=x+p;
                    obj.opt_x=x;
                    if(obj.plots)
                    obj.plot_sum.fvals(end+1)=obj.f.val(x);
                    obj.plot_sum.gnorm(end+1)=norm(g);
                    obj.plot_sum.step_lengths(end+1)=norm(p);
                    end
                end
             end
        end 
        



        


            % algorithm 4.3 page 87
     function lamda=iteratives(~,lamda_0,tr_rad,B,g)

            if(norm(B\(-g))<tr_rad && min(eig(B))>=0)  %check whether the minimizer is within the trust region
                lamda=0;
                return
            end

            lamda_curr=lamda_0;
            lamda_mini=-min(eig(B));%negative smallest eigenvalue of B i.e lamda_l<-lamda1
            D=chol(B+(lamda_curr*eye(2)));
            p=(D.'*D)\(-g);

            while (abs (norm(p) - tr_rad)> 0.001)%checking if norm p is close enough to the TR radius eq 4.37

            D = chol(B + (lamda_curr*eye(2)));  %factorize B

            p=(D.'*D)\(-g);%solve for p

            q=(D.')\p;%find q

            p_norm=norm(p);
            q_norm=norm(q);
            lamda_curr=lamda_curr+((p_norm/q_norm)^2*((p_norm-tr_rad)/tr_rad)); %eq 4.44

            if(lamda_curr<lamda_mini)
                lamda=lamda_curr; 
                return
            end
            end
        lamda=lamda_curr;
            return
        end



function p=piecewise_dogleg(obj,pb,pu,delk)
            if(norm(pb)<delk)
                pk=pb;
            else
                if(norm(pb)>=delk)
                    pk=delk*pu/norm(pu);
                else
                    t=intersection_dog(obj,pb,pu,delk);
                    pk=pu+(t*(pb-pu));
                end
            end
            p=pk;
            return
end

       function t=intersection_dog(~,pb,pu,delk)
%page 75
          % Solving for t from ||pu+(t-1)(pb-pu)^2||^2=del^2
            t=-(2*(pu.')*(pb-pu))+sqrt((2*(pu.')*(pb-pu))^2-(4*(norm(pb-pu)^2)*(norm(pu)^2 -(delk^2))));
            t=t/(2*(norm(pb-pu)^2));
            return 
        end


 function Quasinewton(obj, x_0, submethod)
            x = x_0;                  
            H = eye(2); % approx
            fg = obj.f.grad(x);
            a = 1;
            p = [2,2].';
            
            while (criteria_stop(obj, fg, a*p))
                 if obj.num_step > obj.max_iter
                 break;
                 end
                fg = obj.f.grad(x);
                p = -H * fg;      % eq 6.18
                a_k = linesearch(obj, p, x,"wolfe",true);  % Wolfe LS by algorithm
                
                x_2 = x + a_k*p;
                f2_g = obj.f.grad(x_2);
                
                % Compute s and y
                s = x_2 - x;   %eq 6.5
                y = f2_g - fg;  %eq 6.5
                
                rho = 1 / (y.' * s); %eq 6.14
                
                H = eye(2);
                
                if (submethod == "BFGS")
                   % H = H - (H * (s * s') * H) / (s' * H * s) + (y * y') / (y' * s);
                    H = (eye(2) - rho*s*y.') * H * (eye(2) - rho*y*s.') + (rho * (s*s.')); %eq 6.17 pg 140 %rank 2 form
                elseif (submethod == "DFP")
                   % H=(eye(2) - rho*y*s.') * H * (eye(2) - rho*s*y.') + (rho * (y*y.'));
                    H = H - ((H*(y*y.')*H)/(y.'*H*y)) + ((s*s.')/(y.'*s)); %eq 6.15 pg 139
                end
  
                if (obj.plots)
                    obj.plot_sum.fvals(end+1) = obj.f.val(x);
                    obj.plot_sum.gnorm(end+1) = norm(fg);
                    obj.plot_sum.step_lengths(end+1) = norm(s);
                end                
                x = x_2;
                obj.opt_x = x;  
                obj.num_step = obj.num_step + 1;
            end
        end
        







        function inexactnewton(obj, x_0)%inexact newton using Line-CG to find Quasi Newton direction
    x = x_0;              
    B = eye(2);            % Hessian approximation
    fg = obj.f.grad(x);    
    a = 1;                
    p = [2; 2].';         

    while (criteria_stop(obj, fg, a * p))
         if obj.num_step > obj.max_iter
                 break;
         end
        fg = obj.f.grad(x);
        norm_fg = norm(fg);
        p = LS_CG(fg, B); 
        a = linesearch(obj, p, x,"wolfe",true);
        x_2 = x + a * p;
        s = x_2 - x;
        y = obj.f.grad(x_2) - fg;
       
        % Check for skip/damp condition
        if (y' * s <= 1e-12 * (norm(s) * norm(y)))
        else
            
           B = B - (B * (s * s') * B) / (s' * B * s) + (y * y') / (y' * s); %eq 6.19 pg 140 % BFGS-rank 1 update
          
        end

        x = x_2;
        obj.opt_x = x;

        obj.num_step = obj.num_step + 1;
        if (obj.plots)
            obj.plot_sum.fvals(end + 1) = obj.f.val(x);
            obj.plot_sum.gnorm(end + 1) = norm(fg);
            obj.plot_sum.step_lengths(end + 1) = norm(a * p);
        end
    end


   

    function p = LS_CG(g, H)     
        eta = for_seq(norm_fg,1);
        tol = eta * norm_fg;
         function eta = for_seq(norm_fg, force_type)
   
    switch force_type
        case 1
            eta = min(0.5, sqrt(norm_fg));  %page 169
        case 2
            eta = min(0.25,0.1*sqrt(norm_fg));% random to keep superlinear convergence
        
        otherwise
            error('Invalid force_type');
    end
end
        r = g;       % Residual
        d = -r;      % Search direction 
        z = [0; 0]; 
        while true
            % Check for positive-definiteness of search direction
            if (d.' * H * d <= 0)
                if jth == true
                    p = -g;  
                    return
                else
                    p = z;  
                    return
                end
            end
           
            a = (r.' * r) / (d.' * H * d);  

            z_2 = z + a * d;
            r_2 = r + a * H * d;

            
            if (norm(r_2) < tol)
                p = z_2;
                return
            end

            beta = (r_2.' * r_2) / (r.' * r);
            d = -r_2 + beta * d;
            z = z_2;
            r = r_2;
           % iter_num = false;  
        end
    end
end


      
        
function step_length = linesearch(obj, p, x, method, fallback)
    if method == "back_wol"
        [backtracking_successful, a] = bcktrk(p, x);
        if backtracking_successful
            step_length = a;
            return;
        elseif fallback
            step_length = wolfe_linesearch(p, x);
            return;
        else
            error("Backtracking failed and no fallback specified.");
        end
    elseif method == "wolfe"
        step_length = wolfe_linesearch(p, x);
    else
        error("Unknown line search method specified.");
    end

    function [success, a] = bcktrk(p_k, x_k) %algorithm 3.1

        c1 = 0.02; 
        a_max=2;
        rho = 0.7;       
        a = a_max;    

        imax = 20; %maximum backtrack iterations
        success = false;

        for i = 1:imax
            
            if obj.f.val(x_k + a * p_k) <= obj.f.val(x_k) + c1 * a * dot(obj.f.grad(x_k), p_k) %armijo
                success = true;
                return;
            end
            a = rho * a;
        end
    end


%algoritm 3.5
    function step_length = wolfe_linesearch(p, x) 
        c1 = 0.02;      
        c2 = 0.45;
        a_max =2; 

        phi = @(a) obj.f.val(x + a * p);            
        dphi = @(a) dot(obj.f.grad(x + a * p), p);  

        a_prev = 0;
        a_curr = a_max * 0.5; 
        curr=1;
        while true
            if (phi(a_curr) > phi(0) + c1 * a_curr * dphi(0)) || (phi(a_curr) >= phi(a_prev) && curr>1)                
                step_length = zoom(a_prev, a_curr);
                return;
            end

            if abs(dphi(a_curr)) <= -c2 * dphi(0)
                step_length = a_curr;
                return;
            end

            if dphi(a_curr) >= 0
                step_length = zoom(a_curr, a_prev);
                return;
            end

            a_prev = a_curr;
            a_curr = min(a_curr * 1.2, a_max); % Increase step size, but not beyond a_max
            curr=curr+1;
        end



%algoritm 3.6
        function z = zoom(a_lo, a_hi)

    %         f_lo  = obj.f(x + a_lo*p);
    % g_lo  = dot(obj.f.grad(x + a_lo * p), p);   
    % 
    % f_hi  = obj.f(x + a_hi*p);
    % 
    % a_t = a_lo - (g_lo * (a_hi - a_lo)^2) / (2 * (f_hi - f_lo - g_lo * (a_hi - a_lo)));

            while true
                a_t = (a_lo + a_hi) / 2;%interpolate by bisection
                if abs(a_hi - a_lo) < 1e-16
                    z = a_t; % Interval is too small
                    return;
                end
                phi_a_t = phi(a_t);

                % Check Armijo condition
                if (phi_a_t > phi(0) + c1 * a_t * dphi(0)) || (phi_a_t >= phi(a_lo))
                    a_hi = a_t;
                else
                    dphi_a_t = dphi(a_t);

                    if abs(dphi_a_t) <= -c2 * dphi(0)
                        z = a_t;
                        return;
                    end

                    if dphi_a_t * (a_hi - a_lo) >= 0
                        a_hi = a_lo;
                    end
                    a_lo = a_t;
                end
            end
        end
    end
end

        
        
        
        function stops = criteria_stop(~, fg, step)
            
            %   Stop if gradient is small enough or step is too small.
          stops = (norm(fg) >= 1e-6) && (norm(step) >= 1e-6);

            return
        end
             
        
    end
end
