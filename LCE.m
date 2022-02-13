classdef LCE < handle
    %The class for constitutive model for liquid crystal elastomers, 2D
    
    properties %(Access=protected)
        % Material Parameters:
        MaterialProperties
        % (Internal) Variables in current and previous step
        CurrentState     %Structure including F,l3, Fv,l3v, theta and In, Ine, S for backward Euler
        PreviousState    %Structure including F,l3, Fv,l3v, theta and In, Ine, S for backward Euler
        %State Variables:
        num_Spec  %Number of processes in the spectrum
        p         %Lagrangian Multiplier
        p_eq      %Lagrangian Multiplier, EQ part
        S         %Cauchy Stress
        S_eq      %Cauchy Stress, EQ part
        Fe        %2D elastic deformation gradient
        l3e       %Elastic out of plane stretch
        l0        %Initial normalized step length tensor
        l         %Normalized step length tensor
        l_inv      %Inverse of normalized step length tensor
        w          %Angular Velocity corresponding to Spin Gradient
        %Intermediate Variables
        l_star_E   %Fel0FeT
        l_star     %Fl0FT
        dt         %Time Step  
        %Energy Balance
        D_Dir
        D_Net
        Psi
        %Unused variables for debugging purposes
        debug1 
        debug2 
    end
    
    methods
        function obj = LCE(MaterialProperties,theta)
            %Constuctor
            
            %State Variables in Current and Previous Step
            obj.MaterialProperties = MaterialProperties;
            obj.CurrentState.F = eye(2);
            obj.CurrentState.l3 = 1;
            obj.CurrentState.theta = theta;
            obj.CurrentState.In = 3;
            obj.num_Spec = length(obj.MaterialProperties.muNEQ);
            if obj.num_Spec>0
                obj.CurrentState.Fv = zeros(2,2,obj.num_Spec);
                for ii=1:obj.num_Spec
                    obj.CurrentState.Fv(:,:,ii) = eye(2);
                end
                obj.CurrentState.l3v = ones(obj.num_Spec,1);
                obj.CurrentState.Ine = 3*ones(obj.num_Spec,1);
            end
            obj.CurrentState.S  = 0;
            
            obj.PreviousState = obj.CurrentState;
         
            %Initial State Variables  
            n0 = [cos(theta);sin(theta)];
            obj.l0 =  (1-MaterialProperties.Q)*eye(2)+3*MaterialProperties.Q*(n0*n0');
            obj.l_inv = 1/(1-MaterialProperties.Q)*( eye(2)-3*MaterialProperties.Q/(1+2*MaterialProperties.Q)* (n0*n0') );
            
            %Other (intermediate) variables
            obj.S=zeros(2);
            obj.l_star = obj.l0;
            obj.l3e = ones(obj.num_Spec,1); 
            if obj.num_Spec>0
                for ii=1:obj.num_Spec
                    obj.Fe(:,:,ii) = eye(2);
                end
            end
        end

        
        function ApplyLoad(obj,F,dt)
            %Load to to F after dt, calculates all the internal variables and state variables  
            obj.CurrentState.F = F; %Set current state deformation tensor (2X2)
            obj.dt = dt; %Set time step
            obj.CurrentState.l3 = 1/det(F); %Set lambda_3 (by incompressibility)
            L  = 1/dt*(obj.CurrentState.F-obj.PreviousState.F)/obj.CurrentState.F; %Compute velocity gradient
            OMEGA = 0.5*(L-L'); %Calculate the spin tensor
            obj.w = OMEGA(2,1); %Covenrt to angular velocity
            obj.l_star = F*obj.l0*F';
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Calculate Internal Variables
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if obj.num_Spec>0
                %Compute viscous deformation tensor, the function calls
                %compute_theta inside
                obj.Compute_Fv; 
            else
                %If there is no viscous network deformation branches, 
                %Just compute director angle
                obj.Compute_theta; 
            end
            %%%%%%%%%%%%%%%%%%%%%%%%
            %Calculate Cauchy Stress
            %%%%%%%%%%%%%%%%%%%%%%%%
            %Compute_Stress will not change the object properties and 
            %only return values.
            [Stress, P, Stress_EQ, P_EQ]    =  obj.Compute_Stress;
            %Save the returned values to the object properties
            obj.p_eq=P_EQ;
            obj.S_eq=Stress_EQ-P_EQ*eye(2);
            obj.p = P;
            obj.S = Stress - P*eye(2);
            obj.CurrentState.S = sqrt(sum( sum(obj.S.*obj.S)));
            
            obj.CurrentState.In = sum(sum(  (F*obj.l0*F').*obj.l_inv))+obj.CurrentState.l3^2;
            for ii=1:obj.num_Spec
                obj.CurrentState.Ine(ii) = sum(sum(obj.l_star_E(:,:,ii).*obj.l_inv))+obj.l3e(ii)^2;
            end
            
            %Calculate the terms for the energy balance
            etaD = obj.MaterialProperties.eta_D*exp(-obj.MaterialProperties.ks*obj.PreviousState.S);
            muEQ = obj.MaterialProperties.muEQ;
            muNEQ = obj.MaterialProperties.muNEQ;
            etaN = obj.MaterialProperties.eta_N;
            JmEQ    = obj.MaterialProperties.JmEQ;
            JmNEQ    = obj.MaterialProperties.JmNEQ;
            obj.D_Net = 0;
            obj.Psi = -muEQ/2*JmEQ*log(1-(obj.CurrentState.In-3)/JmEQ);
            obj.D_Dir = etaD * ((obj.CurrentState.theta-obj.PreviousState.theta)/obj.dt-obj.w)^2;
            for ii=1:obj.num_Spec
                Lv = (obj.CurrentState.Fv(:,:,ii)-obj.PreviousState.Fv(:,:,ii))/(obj.CurrentState.Fv(:,:,ii))/obj.dt;
                Lv33 = (obj.CurrentState.l3v(ii)- obj.PreviousState.l3v(ii))/obj.CurrentState.l3v(ii)/obj.dt;
                obj.D_Net = obj.D_Net + etaN*(  sum(sum(Lv.*Lv))+ Lv33^2);
                obj.Psi = obj.Psi-muNEQ(ii)/2*JmNEQ*log(1-(obj.CurrentState.Ine(ii)-3)/JmNEQ)-muNEQ(ii)* log( det(obj.Fe(:,:,ii)) *obj.l3e(ii) );
            end
            
            
        end
        
        
        function [S,P,SEQ,PEQ] = Compute_Stress(obj)
            %Make local copies of material properties and directors
            %so that the expressions is easier to look too bulky
            theta = obj.CurrentState.theta;
            muEQ = obj.MaterialProperties.muEQ;
            muNEQ = obj.MaterialProperties.muNEQ;
            Q = obj.MaterialProperties.Q;
            JmEQ = obj.MaterialProperties.JmEQ;
            JmNEQ    = obj.MaterialProperties.JmNEQ;
            n = [cos(theta);sin(theta)];
            In = obj.PreviousState.In; %Value for Previous State is used!
            
            l_starn = obj.l_star*n;
            k = muEQ/(1-Q)*JmEQ/(JmEQ-In+3);
            c = 1.5*Q/(1+2*Q);
            %Equilibrium Stress
            S = k* ( obj.l_star -c*(n* l_starn' + l_starn*n') );
            P = muEQ*JmEQ/(JmEQ-In+3)* (obj.CurrentState.l3)^2 ; 
            SEQ=S;
            PEQ=P;
            %Nonequilibrium Stress
            for ii =1:obj.num_Spec
                l_star_Eii  = obj.l_star_E(:,:,ii);
                l_star_Enii = l_star_Eii*n; 
                k = muNEQ(ii)/(1-Q); 
                if strcmp(obj.MaterialProperties.type{ii}, 'NC')
                    NGFactor = 1;
                elseif strcmp(obj.MaterialProperties.type{ii}, 'NG')
                    NGFactor = JmNEQ/(JmNEQ-obj.PreviousState.Ine(ii)+3);
                end
                S = S + NGFactor*k*(l_star_Eii - c*(n* l_star_Enii' +  l_star_Enii*n'))-muNEQ(ii)*eye(2);
                P = P + muNEQ(ii)*( NGFactor*obj.l3e(ii)^2 - 1);
            end
        end
        
        
        %Compute_Fv, along with ResFv, computes viscous deformation
        %tensors. Inside the compute_Fv, compute_theta is called to
        %calculate the director angle in each iteration.
        function Compute_Fv(obj)
            JmNEQ = obj.MaterialProperties.JmNEQ;
            %Calculate the in plane part of Fv (2X2)
            fun=@(x)ResFv(obj, x);
            [obj.CurrentState.Fv,~,flag,~] = fsolve(fun, obj.PreviousState.Fv);
            if flag<1
               error("fsolve error when solve 2D Fv");
            end
            
            %Calculate the l3v
            etaN = obj.MaterialProperties.eta_N;
            muNEQ = obj.MaterialProperties.muNEQ;
            l3 = obj.CurrentState.l3;
            for ii=1:obj.num_Spec
                l3vn = obj.PreviousState.l3v(ii);
                if strcmp(obj.MaterialProperties.type{ii}, 'NG')
                    NGFactor =  JmNEQ/(JmNEQ-obj.PreviousState.Ine(ii)+3);
                elseif strcmp(obj.MaterialProperties.type{ii}, 'NC')
                    NGFactor = 1;
                end
                fun = @(x) x - l3vn - muNEQ(ii)/etaN(ii) * (NGFactor * l3^2/x^2 - 1)*x*obj.dt;
                obj.CurrentState.l3v(ii) = fzero(fun, l3vn);
                obj.l3e(ii) = obj.CurrentState.l3/obj.CurrentState.l3v(ii);
            end
        end
        function Res = ResFv(obj,trial)
            JmNEQ =obj.MaterialProperties.JmNEQ;
            for ii=1:obj.num_Spec
                obj.Fe(:,:,ii) = obj.CurrentState.F/trial(:,:,ii);
                obj.l_star_E(:,:,ii) = obj.Fe(:,:,ii)*obj.l0*obj.Fe(:,:,ii)';
            end
            Compute_theta(obj); %CALCULATE DIRECTOR BASED ON GUESS OF FV

            etaN = obj.MaterialProperties.eta_N;
            muNEQ = obj.MaterialProperties.muNEQ;
            
            Res = trial-obj.PreviousState.Fv;
            for ii=1:obj.num_Spec
                if strcmp(obj.MaterialProperties.type{ii}, 'NC')
                    NGFactor = 1;
                elseif strcmp(obj.MaterialProperties.type{ii}, 'NG')
                    NGFactor = JmNEQ/(JmNEQ-obj.PreviousState.Ine(ii)+3);
                end
                Lv = muNEQ(ii)/etaN(ii) * ( NGFactor*(obj.Fe(:,:,ii))'*obj.l_inv*obj.Fe(:,:,ii)*obj.l0 - eye(2) );
                Res(:,:,ii) = Res(:,:,ii) - obj.dt * Lv * trial(:,:,ii);
            end
        end        
        %}
        
        
        function Compute_theta(obj)
            %This can be rewritten using Newton Raphson
            
            %Create local copies of material properties (to make code looks more concise)
            Q = obj.MaterialProperties.Q;
            JmEQ = obj.MaterialProperties.JmEQ;
            JmNEQ    = obj.MaterialProperties.JmNEQ;
            muEQ = obj.MaterialProperties.muEQ;
            muNEQ = obj.MaterialProperties.muNEQ;           
            %Eta is evaluated using stress in the previous time step
            etaD = obj.MaterialProperties.eta_D*exp(-obj.MaterialProperties.ks*obj.PreviousState.S);
            %First neo-classical invariant from previous time step is used.
            IN = obj.PreviousState.In;
            if obj.num_Spec>0
                INe = obj.PreviousState.Ine;
            end
            
            %Solve for theta
            A = JmEQ/(JmEQ-IN+3)*muEQ*(obj.l_star(2,2)-obj.l_star(1,1)); 
            B=  JmEQ/(JmEQ-IN+3)*muEQ*obj.l_star(1,2);
            factor = 1/etaD*3*Q/(1-Q)/(1+2*Q);
            for ii=1:obj.num_Spec
                if strcmp(obj.MaterialProperties.type{ii}, 'NC')
                    NGFactor=1;
                elseif strcmp(obj.MaterialProperties.type{ii}, 'NG')
                    NGFactor = JmNEQ/(JmNEQ-INe(ii)+3);
                end
                A = A + NGFactor *muNEQ(ii)*(obj.l_star_E(2,2,ii)-obj.l_star_E(1,1,ii));
                B = B + NGFactor *muNEQ(ii)*obj.l_star_E(1,2,ii);
            end
            fun = @(x) x-obj.PreviousState.theta - obj.dt*obj.w - factor*obj.dt*(0.5*A*sin(2*x)+B*cos(2*x));
            obj.CurrentState.theta = fzero(fun, obj.PreviousState.theta+rand()*pi/180*5);

            n = [cos(obj.CurrentState.theta); sin(obj.CurrentState.theta)];
            obj.l_inv = 1/(1-Q)*( eye(2)-3*Q/(1+2*Q)* (n*n') );
            obj.l    = (1-Q)*eye(2)+3*Q*(n*n');
        end
        
        
        function UpdateHistory(obj)
            %Update the history variables
            obj.PreviousState = obj.CurrentState;
        end
        
        function ResetHistory(obj)
            %Reset the history variables
            obj.CurrentState = obj.PreviousState;
        end
        
    end
end

