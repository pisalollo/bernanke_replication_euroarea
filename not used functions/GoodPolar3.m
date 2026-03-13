%metodo gnerale per fare sign restricions, non lo useremo


function GR = GoodPolar3(B, normrest, rest, ndraws, horiz)
q = size(B,2);
npar = q*(q-1)/2;
GR = ones(q,q,0);
for j=1:ndraws
    angles = rand(npar,1)*2*pi;
    rotationmat = RotationMatrix(q, angles);
    effects0 = diag(B(abs(normrest),:,2)*rotationmat);
    correctsigns = diag(sign(effects0).*sign(normrest));
    rotationmat = rotationmat*correctsigns;
    for k = horiz
        effects = sign(B(:,:,k)*rotationmat);
        if diag(effects(rest{k}(:,1),abs(rest{k}(:,2))))'*sign(rest{k}(:,2)) < size(rest{k},1)
            k=k-1;
            break
        end
    end
    if k == horiz(end)
        GR = cat(3,GR,rotationmat);
    end
end
