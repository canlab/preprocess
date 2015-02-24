function n = nesttest()
    
    foo = 1:3;

    n.nt1 = @nt1;
    n.nt2 = @nt2;
    
    n = class(n, 'nesttest');
    
    function nt1()
        foo = foo+1;
    end
    
    function nt2()
        disp(foo);
    end
end