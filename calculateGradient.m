%% This script is used to perform calculation on the gradient
function [gradient,err] = calculateGradient(dims,source,trueRec,model,unext,u,uprev)
    recording = zeros(dims.nt,length(dims.recPos),'single');
    gradient  = zeros(dims.ny,dims.nx,'single');
    forwardField = zeros(dims.my,dims.mx,dims.nt,'single'); 
    adjointField = zeros(dims.my,dims.mx,dims.nt,'single');
    err = 0;
    for s = 1:dims.ds:length(dims.srcPos)
        %% Run forward simulation on background model
        uprev(:)=0; u(:)=0; unext(:)=0;
        for t = 1:dims.nt
            % Solve wave equation
            unext = solveWaveEqn(model,dims,dims.srcPos(s),source,t,unext,u,uprev);
            % Update u(x,t)
            uprev = u; u = unext;
            % Check wave equation stability
            r = model*dims.dt/dims.dx + model*dims.dt/dims.dy;
            if r>1
                break
            end
            % Record traces
            recording(t,:) = u(dims.recPos);
            % Save forward field for use in correlation
            forwardField(:,:,t) = u(dims.modely,dims.modelx);
        end
        %% Calculate difference and error
        chi = recording-trueRec(:,:,s);
        % Time reversal
        chi = flipud(chi);
        % Error calculation
        err = err + norm(chi);
        %% Run adjoint simulation
        uprev(:)=0; u(:)=0; unext(:)=0;
        for t = 1:dims.nt
            % Solve wave equation using the difference (chi) as sources
            unext = solveWaveEqn(model,dims,dims.recPos,chi,t,unext,u,uprev);
            % Update u(x,t)
            uprev = u; u = unext;
            % Check wave equation stability
            r = model*dims.dt/dims.dx + model*dims.dt/dims.dy;
            if r>1
                break
            end
            % Save adjoint field for use in correlation
            adjointField(:,:,dims.nt-t+1) = u(dims.modely,dims.modelx);
        end
        %% Correlate
        for t = 2:dims.nt-1
            % Calculate the time derivative of the displacement to gradient.
            dadt=(adjointField(:,:,t+1) - adjointField(:,:,t))/dims.dt;
            dfdt=(forwardField(:,:,t+1) - forwardField(:,:,t-1))/dims.dt;
            gradient(dims.modely,dims.modelx)=gradient(dims.modely,dims.modelx)+dadt.*dfdt;                  
        end
    end
end

