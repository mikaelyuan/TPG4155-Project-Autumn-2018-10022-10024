function [stepLength,newErr] = calculateStepLengthGoldenRatio(dims,model,gradient,source,trueRec)
%% Defining Variables                      
    recordinga = zeros(dims.nt,length(dims.recPos),'single');
    recordingb = zeros(dims.nt,length(dims.recPos),'single');
    recordingc = zeros(dims.nt,length(dims.recPos),'single');
    recordingd = zeros(dims.nt,length(dims.recPos),'single');
    a = 0; 
    b = 256; 
    eps = 10^-2; 
    delta = 10^-2;
    r1 = (sqrt(5)-1)/2; 
    r2 = r1^2; 
    h = b-a;
    c = a+r2*h; 
    d = a+r1*h;
    erra = 0; errb = 0; errc = 0; errd = 0;
    modelnewa = model + a*gradient;
    modelnewb = model + b*gradient;
    modelnewc = model + c*gradient;
    modelnewd = model + d*gradient;
    for s = 1:dims.ds:length(dims.srcPos)
        upreva = zeros(dims.ny,dims.nx,'single');
        ua = zeros(dims.ny,dims.nx,'single');
        unexta = zeros(dims.ny,dims.nx,'single');
        uprevb = zeros(dims.ny,dims.nx,'single');
        ub = zeros(dims.ny,dims.nx,'single');
        unextb = zeros(dims.ny,dims.nx,'single');
        uprevc = zeros(dims.ny,dims.nx,'single');
        uc = zeros(dims.ny,dims.nx,'single');
        unextc = zeros(dims.ny,dims.nx,'single');
        uprevd = zeros(dims.ny,dims.nx,'single');
        ud = zeros(dims.ny,dims.nx,'single');
        unextd = zeros(dims.ny,dims.nx,'single');
        for t = 1:dims.nt
            %  Solve wave equation using test model update
            unexta = solveWaveEqn(modelnewa,dims,dims.srcPos(s),source,t,unexta,ua,upreva);
            unextb = solveWaveEqn(modelnewb,dims,dims.srcPos(s),source,t,unextb,ub,uprevb);
            unextc = solveWaveEqn(modelnewc,dims,dims.srcPos(s),source,t,unextc,uc,uprevc);
            unextd = solveWaveEqn(modelnewd,dims,dims.srcPos(s),source,t,unextd,ud,uprevd);
            upreva = ua; ua = unexta;
            uprevb = ub; ub = unextb;
            uprevc = uc; uc = unextc;
            uprevd = ud; ud = unextd;
            %  Record traces
            recordinga(t,:) = unexta(dims.recPos);
            recordingb(t,:) = unextb(dims.recPos);
            recordingc(t,:) = unextc(dims.recPos);
            recordingd(t,:) = unextd(dims.recPos);
        end
        %% Calculate new error and check against old
        chia = recordinga(:,:)-trueRec(:,:,s);
        chib = recordingb(:,:)-trueRec(:,:,s);
        chic = recordingc(:,:)-trueRec(:,:,s);
        chid = recordingd(:,:)-trueRec(:,:,s);
        erra = erra + norm(chia); 
        errb = errb + norm(chib);
        errc = errc + norm(chic); 
        errd = errd + norm(chid);
    end
    stepLength = b; newErr = errb;
    while (abs(erra-errb)>eps)||(h>delta)
        stepLength = b; newErr = errb;
        if(errc < errd)
            b = d; errb = errd;
            d = c; errd = errc;
            errc = 0;
            h = b-a; c = a+r2*h;
            modelnewc = model + c*gradient;
            for s = 1:dims.ds:length(dims.srcPos)
                uprevc = zeros(dims.ny,dims.nx,'single');
                uc = zeros(dims.ny,dims.nx,'single');
                unextc = zeros(dims.ny,dims.nx,'single');
                for t = 1:dims.nt
                    %  Solve wave equation using test model update
                    unextc = solveWaveEqn(modelnewc,dims,dims.srcPos(s),source,t,unextc,uc,uprevc);
                    uprevc = uc; uc = unextc;
                    %  Record traces
                    recordingc(t,:) = unextc(dims.recPos);
                end
                %% Calculate new error and check against old
                chic = recordingc(:,:)-trueRec(:,:,s);
                errc = errc + norm(chic);
            end
        else
            a = c; 
            erra = errc;
            c = d; 
            errc = errd;
            errd = 0;
            h = b-a; d = a+r1*h;
            modelnewd = model + d*gradient;
            for s = 1:dims.ds:length(dims.srcPos)
                uprevd = zeros(dims.ny,dims.nx,'single');
                ud = zeros(dims.ny,dims.nx,'single');
                unextd = zeros(dims.ny,dims.nx,'single');
                for t = 1:dims.nt
                    %  Solve wave equation using test model update
                    unextd = solveWaveEqn(modelnewd,dims,dims.srcPos(s),source,t,unextd,ud,uprevd);
                    uprevd = ud; ud = unextd;
                    %  Record traces
                    recordingd(t,:) = unextd(dims.recPos);
                end
                %% Calculate new error and check against old
                chid = recordingd(:,:)-trueRec(:,:,s);
                errd = errd + norm(chid);
            end
        end
    end
end