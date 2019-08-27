function  phi2=transport(phi1,u1,N,dx,dt)
    if(size(phi1,1)==N+1 && size(u1,1)==N+1 && length(dx)==1 && length(dt)==1)
        A=zeros(N+1);
        B=zeros(N+1,1);
        A(1,1)=-1;A(1,2)=1;
        for i =2:N
            A(i,i-1)=(-dt/2/dx)*u1(i);
            A(i,i)=1+dt/2/dx*(u1(i+1)-u1(i-1));
            A(i,i+1)=(dt/2/dx)*u1(i);
            B(i)=phi1(i);
        end
        A(N+1,N)=-1;A(N+1,N+1)=1;
        phi2=linsolve(A,B);
    else
        error('Input size is not correct!');
    end
end