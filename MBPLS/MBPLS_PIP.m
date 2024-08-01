function [PIP,tb,t,wb,w] = MBPLS_PIP( xData,yLabel,compNum )

% %   xData: data blocks with B block, xData{b} is the bth data block
% %   yLabel: response data, a column vector
% %   compNumn: number of MB-PLS components
% %   PIP: pathway importance,PIP(b,j) is the eatimated importance of the b
% %        block by the model with j components.
% %   tb:block scores,tb{b}(:,i) is the ith block scores of the bth block  
% %   t: super scores,t(:,i) is the ith super scores of model
% %   wb:block variable weights, wb{b}(:,i) the ith block variable weights of the bth block

Pb = xData;
Y = yLabel;
blockSize = length(Pb);
P = [];
% variable scaling,block scaling
for i = 1:blockSize
    [sampNum,varNum] = size(Pb{i});
    Pb{i} = Pb{i} - ones(sampNum,1)*mean(Pb{i});
    Pb{i} = Pb{i}*diag(1./std(Pb{i}));
    Pb{i} = Pb{i}/sqrt(varNum);
    P = [P Pb{i}];
end
for i = 1:compNum
    V(:,i) = P'*Y;
    T=[];
    for b = 1:blockSize
        wb{b}(:,i) = Pb{b}'*Y;
        wb{b}(:,i) = wb{b}(:,i)/norm(wb{b}(:,i));
        tb{b}(:,i) = Pb{b}*wb{b}(:,i);
        T=[T tb{b}(:,i)];
    end
    w(:,i)=T'*Y;
    w(:,i)=w(:,i)/norm(w(:,i));
    t(:,i)=T*w(:,i);
    for b=1:blockSize
        g{b}(:,i)=Pb{b}'*t(:,i)/(t(:,i)'*t(:,i));
        q(i)=Y'*t(:,i)/(t(:,i)'*t(:,i));
        Pb{b}=Pb{b}-t(:,i)*g{b}(:,i)';
    end
    R=corrcoef(yLabel,t(:,i));
    Y = Y - t(:,i)*q(i)';
    r_square(i,1) = R(2).^2;
    PIP(:,i) = sqrt((blockSize*w.^2*r_square)/sum(r_square));
end

