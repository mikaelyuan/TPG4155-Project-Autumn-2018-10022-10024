%% This script is used to perform calculation on the steplength
function [stepLength,newErr] = calculateStepLength(dims,model,gradient,oldErr,source,trueRec,unext,u,uprev)
%% Depending on how you solve the problem you most probably need more input variables                      
    recording = zeros(dims.nt,length(dims.recPos),'single');
    stepLength = 256;
    newErr = inf;
    while (newErr > oldErr)
        newErr = 0;
        stepLength = stepLength/2;
        %% Test model update calculation based on steplength and gradient
        modelnew = model+stepLength*gradient;
        %% Solve wave equation using test model update
        for s = 1:dims.ds:length(dims.srcPos)
            uprev(:)=0; u(:)=0; unext(:)=0;
            for t = 1:dims.nt
                unext = solveWaveEqn(modelnew,dims,dims.srcPos(s),source,t,unext,u,uprev);
                % Update u(x,t)
                uprev = u; u = unext;
                % Check wave equation stability
                r = model*dims.dt/dims.dx + model*dims.dt/dims.dy;
                if r>1
                    break
                end
                %  Record traces
                recording(t,:) = u(dims.recPos);
            end
            %% Calculate new error and check against old
            chi    = recording-trueRec(:,:,s);
            newErr = newErr+norm(chi);
        end
    end
end