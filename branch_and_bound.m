function branch_and_bound
global x0 flag f0 op n n0;
%线性约束求解器基于内点法或对偶单纯形法，这里为默认选择(对偶单纯形)
%x0为暂时的最优解，flag为标记，op为优化选项，f0为最值，n为初始自变量个数,n0为目标整数在自变量中的位置
%比如 min (3x1 + 2x3) n0应为[1,3]
flag=0;op = optimset('display','off');
%以下为需要手动输入的内容： n,n0,A,b,c;
n = 5;n0 =[1,2];
A =[4,-2,-1,0,0;4,2,0,1,0;0,1,0,0,-1];
b = [1,11,0.5];
c= [-1,-1,0,0,0];
fuc(c,A,b,n);
disp(x0(n0));
disp(f0);
    function [xk,fk]=fuc(c,A,b,n)%更新后的A，b矩阵
    xl = zeros(1,n);
      [xk,fk,exitflag] = linprog(c',[],[],A,b',xl,[],op);
    if exitflag ~= 1 %异常退出
        return
    end
    if flag ~=0 
        if f0<fk
            return;    
        end
    end    %剪枝
    k = 0;
    for i = n0
         if abs(xk(i)-round(xk(i)))>1.0e-4 %判断整数的精度
                break;
         end
         k = k+1;
     end
    if k == size(n0,2)  %整数解
      if flag ==0  %判断x0是否已经有了一个可行整数解
          f0 =fk;
          x0 = xk;
          flag =1;
      elseif f0>fk %判断已知可行解与当前可行解是否更新
            x0 = xk;
            f0 = fk;
      end
    else
        n=n+1;
        A(end+1,end+1) = -1;
        c(end+1) = 0;
        A(end,i) = 1;
        b(end+1) = ceil(xk(i));
        %自变量 n，A,b,的维度扩充；
        fuc(c,A,b,n);% 右迭代
        A(end,end) = 1;
        b(end) = floor(xk(i));
       fuc(c,A,b,n);%  左<=
    end
    
    end


end
