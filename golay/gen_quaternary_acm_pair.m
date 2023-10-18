% to generate quaternary GCM pair of size s1 x s2 constrained by Corollary 2 
% input: m is the length in the first dimension 
%           n is the length in the second dimension
% output: [a, b] is a quaternary GCM pair
% Last modified on April. 30, 2021
% Copyright Communication System Research Laboratory, Fudan University
function [a, b] = gen_quaternary_gcm_pair(m, n)
if m==0 || n==0
    a = []; b = [];
    return;
end
if m==1&& n==1
    a = 1; b = 1; 
    return;
end
em = len_dec(m);
en = len_dec(n);
BBGN = sum(em(1:3)) + sum(en(1:3));
BQGN = sum(em(4:7)) + sum(en(4:7));
assert(BBGN >= BQGN-1, 'Infeasible size');
 
e2 = [1, 1; 1, -1];
e10 =[1 -1 -1 1 -1 1 -1 -1 -1 1;
    1 -1 -1 -1 -1 -1 -1 1 1 -1];
e26 = [-1,1,-1,-1,1,1,-1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,1,1,-1,-1,-1,1,-1,1;
    -1,1,-1,-1,1,1,-1,1,1,1,1,-1,1,-1,1,1,1,1,-1,-1,1,1,1,-1,1,-1];
e3 = [1, 1, -1; 1, 1i, 1];
e5 = [1i, 1i, 1, -1, 1;
    1i, 1, 1, 1i, -1];
e11 = [1, 1i, -1, 1, -1, 1i, -1i, -1, 1i, 1i, 1;
    1, 1, -1i, -1i, -1i, 1, 1, 1i, -1, 1, -1];
e13 = [1, 1, 1, 1i, -1, 1, 1, -1i, 1, -1, 1, -1i, 1i;
    1, 1i, -1, -1, -1, 1i, -1, 1, 1, -1i, -1, 1, -1i];
elements = {e2, e10, e26, e3, e5, e11, e13};
expo = [em; en];
index = 1;
glues = [];
for i = 1:3
    for j = 1:2
        while expo(j, i) > 0
            if j ==1
                glues{index} = {elements{i}(1, :).', elements{i}(2, :).'};
            else
                glues{index} = {elements{i}(1, :), elements{i}(2, :)};
            end
            expo(j, i) = expo(j, i)-1;
            index = index+1;
        end
    end
end

papers = [];
index = 1;
for i = 4:7
    for j = 1:2
        while expo(j, i) > 0
            if j ==1
                papers{index} = {elements{i}(1, :).', elements{i}(2, :).'};
            else
                papers{index} = {elements{i}(1, :), elements{i}(2, :)};
            end
            expo(j, i) = expo(j, i)-1;
            index = index+1;
        end
    end
end

if length(papers) > 1
    a = papers{1}{1};
    b = papers{1}{2};
    for i = 1:length(papers)-1
        c = papers{i+1}{1};
        d = papers{i+1}{2};
        x = glues{i}{1};
        y = glues{i}{2};
        [a, b] = golay_glue(x, y, a, b, c, d);
    end
    for j = 1:length(glues)-i
        x = glues{i+j}{1};
        y = glues{i+j}{2};
        [a, b] = golay_glue(x, y, a, b, 1, 1);
    end
else
    if length(papers) == 1
        a = papers{1}{1};
        b = papers{1}{2};
    else 
        a = 1; b = 1;
    end
    for j = 1:length(glues)
        x = glues{j}{1};
        y = glues{j}{2};
        [a, b] = golay_glue(x, y, a, b, 1, 1);
    end
end
end

function expo = len_dec(num)
init = num;
base = [2, 3, 5, 11, 13];
expo = zeros(1, 5);
for i = 1:length(base)
    [expo(i), num] = expo_dec(num, base(i));
end
if num ~= 1
    error(sprintf('%d is undecomposable', init));
end

expo10 = min(expo(1), expo(3));
expo26 = min(expo(1)-expo10, expo(5));
expo = [expo(1)-expo10-expo26, expo10, expo26, ...
    expo(2), expo(3)-expo10, expo(4), expo(5)-expo26];
end


function [expo, num] = expo_dec(num, base)
expo = 0;
while mod(num, base) == 0
    num = num / base;
    expo = expo + 1;
end
end

function [e, f] = golay_glue(x, y, a, b, c, d)
u = 1/4 * (x+flip_all(x)+y-flip_all(y));
v = 1/4 * (x-flip_all(x)+y+flip_all(y));
s = kron(a, u) + kron(b, v);
t = kron(a, flip_all(conj(v))) - kron(b, flip_all(conj(u)));
e = kron(c, s) + kron(flip_all(conj(d)), t);
f = kron(flip_all(conj(c)), t) - kron(d, s);
end

function a = flip_all(a)
sa = size(a);
for i = 1:length(sa)
    a = flip(a, i);
end
end
