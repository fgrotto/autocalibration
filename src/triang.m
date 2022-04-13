function X = triang(P1, P2, m1, m2) 
    x1 = m1(1,1);
    y1 = m1(2,1);
    x2 = m2(1,1);
    y2 = m2(2,1);
    A = [(P1(1,:) - x1*(P1(3,:)));
        (P1(2,:) - y1*(P1(3,:)));
        (P2(1,:) - x2*(P2(3,:)));
        (P2(2,:) - y2*(P2(3,:)))];

    [~,~,V] = svd(A,0);
    M = V(:,end);

    X = [M(1,:)./M(4,:) M(2,:)./M(4,:) M(3,:)./M(4,:)];
end