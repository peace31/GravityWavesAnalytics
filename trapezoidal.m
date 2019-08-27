function y=trapezoidal(f,N,dx)
    if(size(f,1)==N+1 && length(dx)==1)
        y=0;
        for i=1:N
            y=y+(f(i)+f(i+1))/2;
        end
    else
        error('Input size is not correct!');
    end
end