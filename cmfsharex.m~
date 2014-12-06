function [f1,f2, s1 , g, u, v] = cmfsharex(x1, x2, f1, f2, g, eqnum, iter, test,s,w1,w2)

%beta cmf ||X1 - (F;F1)S1G'|| + ||X2 - (F;F2)(U*S2*V)G'||
%X=FSG, F and G are normalized \sum F_j = 1
%author:antonyxie

if nargin < 10
    w1 = ones(size(x1,1), size(x1,2));
    w2 = ones(size(x2,1), size(x2,2));
end
%s = abs(pinv(f'*f)*f'*x*g*pinv(g'*g));
s1 = rand(size(f1,2), size(g,2));
s2 = s1;
s1 = s;
s2 =s;
u = eye(size(f1,2));
v = eye(size(g,2));

f2(1:eqnum,:) = f1(1:eqnum,:);
f = f1(1:eqnum,:);
fc = f1(eqnum+1:size(f1,1), :);
fd = f2(eqnum+1:size(f2,1), :);

IF1 = eye(size(f1,1));
Icc = IF1(:, 1:eqnum);
Ic = IF1(:, eqnum+1:size(f1,1));

IF2 = eye(size(f2,1));
Icd = IF2(:, 1:eqnum);
Id = IF2(:, eqnum+1:size(f2,1));

%iter = 30;
a = 0.0001;
b = 1;

prec = zeros(10,1);

biccx1 = sparse(b*Icc'*(w1.*x1));
icdx2 = sparse(Icd'*(w2.*x2));
bicx1 = sparse(b*Ic'*(w1.*x1));
bx1 = sparse(b*(w1.*x1)');
idx2 = sparse(Id'*(w2.*x2));
biccicc = sparse(b*Icc'*Icc);
biccic = sparse(b*Icc'*Ic);
bicicc = sparse(b*Ic'*Icc);
bicic = sparse(b*Ic'*Ic);
for j = 1:10
    for i = 1:iter

        %#Ac = Icc*f*s1*g';
        %#Bc = Ic*fc*s1*g';
        %#Ad = Icd*f*s2*g';
        %#Bd = Id*fd*s2*g';
        %#Icc'*x1*s1*g;
        s2 = s1;
        sg = s1*g'*g*s1';
        sg2 = s2*g'*g*s2';
        tic
        
        %( b*Icc'*(Icc*f + Ic*fc)) *s1*g'*g*s1';
        f = f .* sqrt((biccx1*g*s1' + icdx2*g*s2') ./  ( (biccicc*f + biccic*fc)*sg + ( Icd'*(Icd*f + (Id*fd)))*sg2  + a*f)) + 10^-18;
        %f = f ./ (sum(f')' * ones(1, size(f,2)));
        toc
        %#fc
        tic
        fc = fc .* sqrt((bicx1*g*s1') ./ ( (bicicc*f + bicic*fc)*sg + a*fc)) + 10^-18;
        toc
        %fc = fc ./ (sum(fc')' * ones(1, size(fc,2)));
        %%#f1(1:eq
        %#fd
        tic
        fd = fd .* sqrt((idx2*g*s2') ./ ((Id'*(Icd*f + Id*fd)) *sg2 + a*fd)) + 10^-18;
        toc
       % fd = fd ./ (sum(fd')' * ones(1, size(fd,2)));

        %#size(f)
        %fd
       % s1
        f1 = [f;fc];
        %#size(f1)
        f2 = [f;fd];
        %#size(f2)
        %f1 = f1 .* sqrt((x1*g*s1')./(f1*s1*g'*g*s1'));
        %normalize
        %f1 = f1 ./ (sum(f1')' * ones(1, size(f1,2)));

        %f2 = f2 .* sqrt((x2*g*s2')./(f2*s2*g'*g*s2'));
        %normalize
        %f2 = f2 ./ (sum(f2')' * ones(1, size(f2,2)));

        %
        %f2(1:eqnum,:) = f1(1:eqnum,:);
        tic
        %#x1'*f1*s1
        g = g .* sqrt( (bx1*f1*s1 + x2'*f2*s2) ./ ((g*(b*s1'*f1'*f1*s1 + s2'*f2'*f2*s2) + a*g) )) + 10^-18;
       % g = g ./ (sum(g')' * ones(1, size(g,2)));
        toc

            
        tic
        s1 = s1 .* sqrt((b*f1'*x1*g + (f2)'*x2*(g))./((b*f1'*f1*s1 + f2'*f2*s2)*g'*g) ) + 10^-18;
        toc
        s2 = s1;
    %for l = 1:10
        u = u .* sqrt( (f2'*x2*g*v'*s2') ./ (f2'*f2*u*s2*v*g'*g*v'*s2')) ;
       % u = u * (size(u,1) / sum(sum(u)));

        %v = v .* sqrt( (s2'*f2'*x2*g) ./ (s2'*f2'*f2*s2*v*g'*g));
        v = v * (size(v,1) / sum(sum(v)));
    %end

        %#g(1,1)
        %#f'*f
        %#s1(1,1)
        fsg = f1*s1*g';
        fsg2 = f2*(u*s2*v)*g';
        %fsg3 = f2*s2*g';wcmfsharex(trainpv(1:8000,:), train, ff, fff, gt, 2100, 20, ones(8000,1791), ones(2520,1791), test);


        %err1 = sqrt(trace((x1-fsg)'*(x1-fsg)) )
        %err2 =  sqrt(trace((x2-fsg2)'*(x2-fsg2)))
        %err3 = 	sqrt(trace((x2-fsg3)'*(x2-fsg3)))
        i
        %err = sqrt(trace((x1-fsg)'*(x1-fsg)) + trace((x2-fsg2)'*(x2-fsg2)))
        %errr = sqrt(b*sum(sum((x1-fsg).^2)) + sum(sum((x2-fsg2).^2)) + a*(sum(sum(f1.^2)) + sum(sum(f2.^2)) + sum(sum(g.^2)) ))

    end
    er = testprec( fsg2 - x2*1000, test);
    prec(j) = er;
    if  j > 1 & prec(j) < prec(j-1) 
        j
        %break
    end
end
