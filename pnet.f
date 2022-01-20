      subroutine plognet (parm,no,ni,nc,x,y,g,jd,vp,cl,ne,nx,nlam,flmin,u   1533 
     *lam,thr,isd,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,nc
     *p,pp, jerr)
      implicit double precision(a-h,o-z)                                   1534
      double precision x(no,ni),y(ncp,max(2,nc)),g(ncp,nc),vp(ni),         1535 
     *ulam(nlam)
      double precision ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl   1536 
     *(2,ni)
      integer jd(*),ia(nx),nin(nlam)                                       1537
      double precision, dimension (:), allocatable :: xm,xs,ww,vq,xv            
      integer, dimension (:), allocatable :: ju
      
      integer intr
      double precision x2(ncp,ni)
      integer pp(ncp,2)
      intr = 0 

c     print *, "flmin = ", flmin
      x2 =  x(pp(:,1),:) - x(pp(:,2),:)

      if(maxval(vp) .gt. 0.0)goto 12221                                    1541
      jerr=10000                                                           1541
      return                                                               1541
12221 continue                                                             1542
      allocate(ww(1:ncp),stat=jerr)                                        1543
      if(jerr.ne.0) return                                                 1544
      allocate(ju(1:ni),stat=jerr)                                         1545
      if(jerr.ne.0) return                                                 1546
      allocate(vq(1:ni),stat=jerr)                                         1547
      if(jerr.ne.0) return                                                 1548
      allocate(xm(1:ni),stat=jerr)                                         1549
      if(jerr.ne.0) return                                                 1550
      if(kopt .ne. 2)goto 12241                                            1550
      allocate(xv(1:ni),stat=jerr)                                         1550
      if(jerr.ne.0) return                                                 1550
12241 continue                                                             1551
      if(isd .le. 0)goto 12261                                             1551
      allocate(xs(1:ni),stat=jerr)                                         1551
      if(jerr.ne.0) return                                                 1551
12261 continue                                                             1553
      call chkvars(ncp,ni,x2,ju)                                          1554
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 1555
      if(maxval(ju) .gt. 0)goto 12281                                      1555
      jerr=7777                                                            1555
      return                                                               1555
12281 continue                                                             1556
      vq=max(0d0,vp)                                                       1556
      vq=vq*ni/sum(vq)                                                     1557
12290 do 12291 i=1,ncp                                                     1557
      ww(i)=sum(y(i,:))                                                    1557
      if(ww(i).gt.0.0) y(i,:)=y(i,:)/ww(i)                                 1557
12291 continue                                                             1558
12292 continue                                                             1558
      sw=sum(ww)                                                           1558
      ww=ww/sw                                                             1559
      if(nc .ne. 1)goto 12311                                              1559
      call lstandard1(ncp,ni,x2,ww,ju,isd,intr,xm,xs)                     1560
      if(isd .le. 0)goto 12331                                             1560
12340 do 12341 j=1,ni                                                      1560
      cl(:,j)=cl(:,j)*xs(j)                                                1560
12341 continue                                                             1560
12342 continue                                                             1560
12331 continue                                                             1561
      call lognet2n(parm,ncp,ni,x2,y(:,1),g(:,1),ww,ju,vq,cl,ne,nx,nlam, 1563 
     *flmin,ulam,  thr,isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,
     *alm,nlp,jerr)
      goto 12301                                                           1564
12311 if(kopt .ne. 2)goto 12351                                            1564
c     call multlstandard1(ncp,ni,x2,ww,ju,isd,intr,xm,xs,xv)              1565
c     if(isd .le. 0)goto 12371                                             1565
12380 do 12381 j=1,ni                                                      1565
      cl(:,j)=cl(:,j)*xs(j)                                                1565
12381 continue                                                             1565
12382 continue                                                             1565
12371 continue                                                             1566
c     call multlognetn(parm,ncp,ni,nc,x2,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,1568 
c    *ulam,thr,  intr,maxit,xv,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
c     goto 12391                                                           1569
12351 continue                                                             1569
      call lstandard1(ncp,ni,x2,ww,ju,isd,intr,xm,xs)                     1570
      if(isd .le. 0)goto 12411                                             1570
12420 do 12421 j=1,ni                                                      1570
      cl(:,j)=cl(:,j)*xs(j)                                                1570
12421 continue                                                             1570
12422 continue                                                             1570
12411 continue                                                             1571
      call lognetn(parm,ncp,ni,nc,x2,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,     1573 
     *ulam,thr,isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,
     *jerr)
12391 continue                                                             1574
12301 continue                                                             1574
      if(jerr.gt.0) return                                                 1574
      dev0=2.0*sw*dev0                                                     1575
12430 do 12431 k=1,lmu                                                     1575
      nk=nin(k)                                                            1576
12440 do 12441 ic=1,nc                                                     1576
      if(isd .le. 0)goto 12461                                             1576
12470 do 12471 l=1,nk                                                      1576
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      1576
12471 continue                                                             1576
12472 continue                                                             1576
12461 continue                                                             1577
      if(intr .ne. 0)goto 12491                                            1577
      a0(ic,k)=0.0                                                         1577
      goto 12501                                                           1578
12491 continue                                                             1578
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            1578
12501 continue                                                             1579
12481 continue                                                             1579
12441 continue                                                             1580
12442 continue                                                             1580
12431 continue                                                             1581
12432 continue                                                             1581
      deallocate(ww,ju,vq,xm)                                              1581
      if(isd.gt.0) deallocate(xs)                                          1582
      if(kopt.eq.2) deallocate(xv)                                         1583
      return                                                               1584
      end                                                                  1585
      
      subroutine chkvars(no,ni,x,ju)                                       1130
      implicit double precision(a-h,o-z)                                   1131
      double precision x(no,ni)                                            1131
      integer ju(ni)                                                       1132
11060 do 11061 j=1,ni                                                      1132
      ju(j)=0                                                              1132
      t=x(1,j)                                                             1133
11070 do 11071 i=2,no                                                      1133
      if(x(i,j).eq.t)goto 11071                                            1133
      ju(j)=1                                                              1133
      goto 11072                                                           1133
11071 continue                                                             1134
11072 continue                                                             1134
11061 continue                                                             1135
11062 continue                                                             1135
      return                                                               1136
      end                                                                  1137
      
      subroutine lstandard1 (no,ni,x,w,ju,isd,intr,xm,xs)                  1586
      implicit double precision(a-h,o-z)                                   1587
      double precision x(no,ni),w(no),xm(ni),xs(ni)                        1587
      integer ju(ni)                                                       1588
      if(intr .ne. 0)goto 12521                                            1589
12530 do 12531 j=1,ni                                                      1589
      if(ju(j).eq.0)goto 12531                                             1589
      xm(j)=0.0                                                            1590
      if(isd .eq. 0)goto 12551                                             1590
      vc=dot_product(w,x(:,j)**2)-dot_product(w,x(:,j))**2                 1591
      xs(j)=sqrt(vc)                                                       1591
      x(:,j)=x(:,j)/xs(j)                                                  1592
12551 continue                                                             1593
12531 continue                                                             1594
12532 continue                                                             1594
      return                                                               1595
12521 continue                                                             1596
12560 do 12561 j=1,ni                                                      1596
      if(ju(j).eq.0)goto 12561                                             1597
      xm(j)=dot_product(w,x(:,j))                                          1597
      x(:,j)=x(:,j)-xm(j)                                                  1598
      if(isd .le. 0)goto 12581                                             1598
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                                 1598
      x(:,j)=x(:,j)/xs(j)                                                  1598
12581 continue                                                             1599
12561 continue                                                             1600
12562 continue                                                             1600
      return                                                               1601
      end                                                                  1602
      subroutine lognet2n(parm,no,ni,x,y,g,w,ju,vp,cl,ne,nx,nlam,flmin,u   1623 
     *lam,shri,  isd,intr,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jer
     *r)
      implicit double precision(a-h,o-z)                                   1624
      double precision x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam),cl(2   1625 
     *,ni)
      double precision a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)             1626
      integer ju(ni),m(nx),kin(nlam)                                       1627
      double precision, dimension (:), allocatable :: b,bs,v,r,xv,q,ga          
      integer, dimension (:), allocatable :: mm,ixx                             
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               1632
      allocate(b(0:ni),stat=jerr)                                          1633
      if(jerr.ne.0) return                                                 1634
      allocate(xv(1:ni),stat=jerr)                                         1635
      if(jerr.ne.0) return                                                 1636
      allocate(ga(1:ni),stat=jerr)                                         1637
      if(jerr.ne.0) return                                                 1638
      allocate(bs(0:ni),stat=jerr)                                         1639
      if(jerr.ne.0) return                                                 1640
      allocate(mm(1:ni),stat=jerr)                                         1641
      if(jerr.ne.0) return                                                 1642
      allocate(ixx(1:ni),stat=jerr)                                        1643
      if(jerr.ne.0) return                                                 1644
      allocate(r(1:no),stat=jerr)                                          1645
      if(jerr.ne.0) return                                                 1646
      allocate(v(1:no),stat=jerr)                                          1647
      if(jerr.ne.0) return                                                 1648
      allocate(q(1:no),stat=jerr)                                          1649
      if(jerr.ne.0) return                                                 1650
      fmax=log(1.0/pmin-1.0)                                               1650
      fmin=-fmax                                                           1650
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      1651
      bta=parm                                                             1651
      omb=1.0-bta                                                          1652
      q0=dot_product(w,y)                                                  1652
      if(q0 .gt. pmin)goto 12681                                           1652
      jerr=8001                                                            1652
      return                                                               1652
12681 continue                                                             1653
      if(q0 .lt. 1.0-pmin)goto 12701                                       1653
      jerr=9001                                                            1653
      return                                                               1653
12701 continue                                                             1654
      if(intr.eq.0.0) q0=0.5                                               1655
      ixx=0                                                                1655
      al=0.0                                                               1655
      bz=0.0                                                               1655
      if(intr.ne.0) bz=log(q0/(1.0-q0))                                    1656
      if(nonzero(no,g) .ne. 0)goto 12721                                   1656
      vi=q0*(1.0-q0)                                                       1656
      b(0)=bz                                                              1656
      v=vi*w                                                               1657
      r=w*(y-q0)                                                           1657
      q=q0                                                                 1657
      xmz=vi                                                               1657
      dev1=-(bz*q0+log(1.0-q0))                                            1658
      goto 12731                                                           1659
12721 continue                                                             1659
      b(0)=0.0                                                             1660
      if(intr .eq. 0)goto 12751                                            1660
      b(0)=azero(no,y,g,w,jerr)                                            1660
      if(jerr.ne.0) return                                                 1660
12751 continue                                                             1661
      q=1.0/(1.0+exp(-b(0)-g))                                             1661
      v=w*q*(1.0-q)                                                        1661
      r=w*(y-q)                                                            1661
      xmz=sum(v)                                                           1662
      dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)))                        1663
12731 continue                                                             1664
12711 continue                                                             1664
      if(kopt .le. 0)goto 12771                                            1665
      if(isd .le. 0 .or. intr .eq. 0)goto 12791                            1665
      xv=0.25                                                              1665
      goto 12801                                                           1666
12791 continue                                                             1666
12810 do 12811 j=1,ni                                                      1666
      if(ju(j).ne.0) xv(j)=0.25*dot_product(w,x(:,j)**2)                   1666
12811 continue                                                             1666
12812 continue                                                             1666
12801 continue                                                             1667
12781 continue                                                             1667
12771 continue                                                             1668
      dev0=dev1                                                            1669
12820 do 12821 i=1,no                                                      1669
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        1670
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              1671
12821 continue                                                             1673
12822 continue                                                             1673
      alf=1.0                                                              1675
      if(flmin .ge. 1.0)goto 12841                                         1675
      eqs=max(eps,flmin)                                                   1675
      alf=eqs**(1.0/(nlam-1))                                              1675
12841 continue                                                             1676
      m=0                                                                  1676
      mm=0                                                                 1676
      nlp=0                                                                1676
      nin=nlp                                                              1676
      mnl=min(mnlam,nlam)                                                  1676
      bs=0.0                                                               1676
      b(1:ni)=0.0                                                          1677
      shr=shri*dev0                                                        1678
12850 do 12851 j=1,ni                                                      1678
      if(ju(j).eq.0)goto 12851                                             1678
      ga(j)=abs(dot_product(r,x(:,j)))                                     1678
12851 continue                                                             1679
12852 continue                                                             1679
12860 do 12861 ilm=1,nlam                                                  1679
      al0=al                                                               1680
      if(flmin .lt. 1.0)goto 12881                                         1680
      al=ulam(ilm)                                                         1680
      goto 12871                                                           1681
12881 if(ilm .le. 2)goto 12891                                             1681
      al=al*alf                                                            1681
      goto 12871                                                           1682
12891 if(ilm .ne. 1)goto 12901                                             1682
      al=big                                                               1682
      goto 12911                                                           1683
12901 continue                                                             1683
      al0=0.0                                                              1684
12920 do 12921 j=1,ni                                                      1684
      if(ju(j).eq.0)goto 12921                                             1684
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1684
12921 continue                                                             1685
12922 continue                                                             1685
      al0=al0/max(bta,1.0d-3)                                              1685
      al=alf*al0                                                           1686
12911 continue                                                             1687
12871 continue                                                             1687
      al2=al*omb                                                           1687
      al1=al*bta                                                           1687
      tlam=bta*(2.0*al-al0)                                                1688
12930 do 12931 k=1,ni                                                      1688
      if(ixx(k).eq.1)goto 12931                                            1688
      if(ju(k).eq.0)goto 12931                                             1689
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     1690
12931 continue                                                             1691
12932 continue                                                             1691
10880 continue                                                             1692
12940 continue                                                             1692
12941 continue                                                             1692
      bs(0)=b(0)                                                           1692
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                1693
      if(kopt .ne. 0)goto 12961                                            1694
12970 do 12971 j=1,ni                                                      1694
      if(ixx(j).gt.0) xv(j)=dot_product(v,x(:,j)**2)                       1694
12971 continue                                                             1695
12972 continue                                                             1695
12961 continue                                                             1696
12980 continue                                                             1696
12981 continue                                                             1696
      nlp=nlp+1                                                            1696
      dlx=0.0                                                              1697
12990 do 12991 k=1,ni                                                      1697
      if(ixx(k).eq.0)goto 12991                                            1698
      bk=b(k)                                                              1698
      gk=dot_product(r,x(:,k))                                             1699
      u=gk+xv(k)*b(k)                                                      1699
      au=abs(u)-vp(k)*al1                                                  1700
      if(au .gt. 0.0)goto 13011                                            1700
      b(k)=0.0                                                             1700
      goto 13021                                                           1701
13011 continue                                                             1702
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))          1703
13021 continue                                                             1704
13001 continue                                                             1704
      d=b(k)-bk                                                            1704
      if(abs(d).le.0.0)goto 12991                                          1704
      dlx=max(dlx,xv(k)*d**2)                                              1705
      r=r-d*v*x(:,k)                                                       1706
      if(mm(k) .ne. 0)goto 13041                                           1706
      nin=nin+1                                                            1706
      if(nin.gt.nx)goto 12992                                              1707
      mm(k)=nin                                                            1707
      m(nin)=k                                                             1708
13041 continue                                                             1709
12991 continue                                                             1710
12992 continue                                                             1710
      if(nin.gt.nx)goto 12982                                              1711
      d=0.0                                                                1711
      if(intr.ne.0) d=sum(r)/xmz                                           1712
      if(d .eq. 0.0)goto 13061                                             1712
      b(0)=b(0)+d                                                          1712
      dlx=max(dlx,xmz*d**2)                                                1712
      r=r-d*v                                                              1712
13061 continue                                                             1713
      if(dlx.lt.shr)goto 12982                                             1713
      if(nlp .le. maxit)goto 13081                                         1713
      jerr=-ilm                                                            1713
      return                                                               1713
13081 continue                                                             1714
13090 continue                                                             1714
13091 continue                                                             1714
      nlp=nlp+1                                                            1714
      dlx=0.0                                                              1715
13100 do 13101 l=1,nin                                                     1715
      k=m(l)                                                               1715
      bk=b(k)                                                              1716
      gk=dot_product(r,x(:,k))                                             1717
      u=gk+xv(k)*b(k)                                                      1717
      au=abs(u)-vp(k)*al1                                                  1718
      if(au .gt. 0.0)goto 13121                                            1718
      b(k)=0.0                                                             1718
      goto 13131                                                           1719
13121 continue                                                             1720
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))          1721
13131 continue                                                             1722
13111 continue                                                             1722
      d=b(k)-bk                                                            1722
      if(abs(d).le.0.0)goto 13101                                          1722
      dlx=max(dlx,xv(k)*d**2)                                              1723
      r=r-d*v*x(:,k)                                                       1724
13101 continue                                                             1725
13102 continue                                                             1725
      d=0.0                                                                1725
      if(intr.ne.0) d=sum(r)/xmz                                           1726
      if(d .eq. 0.0)goto 13151                                             1726
      b(0)=b(0)+d                                                          1726
      dlx=max(dlx,xmz*d**2)                                                1726
      r=r-d*v                                                              1726
13151 continue                                                             1727
      if(dlx.lt.shr)goto 13092                                             1727
      if(nlp .le. maxit)goto 13171                                         1727
      jerr=-ilm                                                            1727
      return                                                               1727
13171 continue                                                             1728
      goto 13091                                                           1729
13092 continue                                                             1729
      goto 12981                                                           1730
12982 continue                                                             1730
      if(nin.gt.nx)goto 12942                                              1731
13180 do 13181 i=1,no                                                      1731
      fi=b(0)+g(i)                                                         1732
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin)),x(i,m(1:nin)))            1733
      if(fi .ge. fmin)goto 13201                                           1733
      q(i)=0.0                                                             1733
      goto 13191                                                           1733
13201 if(fi .le. fmax)goto 13211                                           1733
      q(i)=1.0                                                             1733
      goto 13221                                                           1734
13211 continue                                                             1734
      q(i)=1.0/(1.0+exp(-fi))                                              1734
13221 continue                                                             1735
13191 continue                                                             1735
13181 continue                                                             1736
13182 continue                                                             1736
      v=w*q*(1.0-q)                                                        1736
      xmz=sum(v)                                                           1736
      if(xmz.le.vmin)goto 12942                                            1736
      r=w*(y-q)                                                            1737
      if(xmz*(b(0)-bs(0))**2 .ge. shr)goto 13241                           1737
      ix=0                                                                 1738
13250 do 13251 j=1,nin                                                     1738
      k=m(j)                                                               1739
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 13251                           1739
      ix=1                                                                 1739
      goto 13252                                                           1740
13251 continue                                                             1741
13252 continue                                                             1741
      if(ix .ne. 0)goto 13271                                              1742
13280 do 13281 k=1,ni                                                      1742
      if(ixx(k).eq.1)goto 13281                                            1742
      if(ju(k).eq.0)goto 13281                                             1743
      ga(k)=abs(dot_product(r,x(:,k)))                                     1744
      if(ga(k) .le. al1*vp(k))goto 13301                                   1744
      ixx(k)=1                                                             1744
      ix=1                                                                 1744
13301 continue                                                             1745
13281 continue                                                             1746
13282 continue                                                             1746
      if(ix.eq.1) go to 10880                                              1747
      goto 12942                                                           1748
13271 continue                                                             1749
13241 continue                                                             1750
      goto 12941                                                           1751
12942 continue                                                             1751
      if(nin .le. nx)goto 13321                                            1751
      jerr=-10000-ilm                                                      1751
      goto 12862                                                           1751
13321 continue                                                             1752
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                1752
      kin(ilm)=nin                                                         1753
      a0(ilm)=b(0)                                                         1753
      alm(ilm)=al                                                          1753
      lmu=ilm                                                              1754
      devi=dev2(no,w,y,q,pmin)                                             1755
      dev(ilm)=(dev1-devi)/dev0                                            1755
      if(xmz.le.vmin)goto 12862                                            1756
      if(ilm.lt.mnl)goto 12861                                             1756
      if(flmin.ge.1.0)goto 12861                                           1757
      me=0                                                                 1757
13330 do 13331 j=1,nin                                                     1757
      if(a(j,ilm).ne.0.0) me=me+1                                          1757
13331 continue                                                             1757
13332 continue                                                             1757
      if(me.gt.ne)goto 12862                                               1758
      if(dev(ilm).gt.devmax)goto 12862                                     1758
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 12862                             1759
12861 continue                                                             1760
12862 continue                                                             1760
      g=log(q/(1.0-q))                                                     1761
      deallocate(b,bs,v,r,xv,q,mm,ga,ixx)                                  1762
      return                                                               1763
      end                                                                  1764
      subroutine lognetn(parm,no,ni,nc,x,y,g,w,ju,vp,cl,ne,nx,nlam,flmin   1799 
     *,ulam,shri,  isd,intr,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,j
     *err)
      implicit double precision(a-h,o-z)                                   1800
      double precision x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam   1801 
     *)
      double precision a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(   1802 
     *2,ni)
      integer ju(ni),m(nx),kin(nlam)                                       1803
      double precision, dimension (:,:), allocatable :: q                       
      double precision, dimension (:), allocatable :: sxp,sxpl                  
      double precision, dimension (:), allocatable :: di,v,r,ga                 
      double precision, dimension (:,:), allocatable :: b,bs,xv                 
      integer, dimension (:), allocatable :: mm,is,ixx                          
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      if(jerr.ne.0) return                                                      
      allocate(xv(1:ni,1:nc),stat=jerr)                                         
      if(jerr.ne.0) return                                                                             
      allocate(bs(0:ni,1:nc),stat=jerr)                                         
      if(jerr.ne.0) return                                                                             
      allocate(q(1:no,1:nc),stat=jerr)                                          
      if(jerr.ne.0) return                                                              
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               1818
      exmn=-exmx                                                           1819
      allocate(r(1:no),stat=jerr)                                          1820
      if(jerr.ne.0) return                                                 1821
      allocate(v(1:no),stat=jerr)                                          1822
      if(jerr.ne.0) return                                                 1823
      allocate(mm(1:ni),stat=jerr)                                         1824
      if(jerr.ne.0) return                                                 1825
      allocate(is(1:max(nc,ni)),stat=jerr)                                 1826
      if(jerr.ne.0) return                                                 1827
      allocate(sxp(1:no),stat=jerr)                                        1828
      if(jerr.ne.0) return                                                 1829
      allocate(sxpl(1:no),stat=jerr)                                       1830
      if(jerr.ne.0) return                                                 1831
      allocate(di(1:no),stat=jerr)                                         1832
      if(jerr.ne.0) return                                                 1833
      allocate(ga(1:ni),stat=jerr)                                         1834
      if(jerr.ne.0) return                                                 1835
      allocate(ixx(1:ni),stat=jerr)                                        1836
      if(jerr.ne.0) return                                                 1837
      pmax=1.0-pmin                                                        1837
      emin=pmin/pmax                                                       1837
      emax=1.0/emin                                                        1838
      pfm=(1.0+pmin)*pmin                                                  1838
      pfx=(1.0-pmin)*pmax                                                  1838
      vmin=pfm*pmax                                                        1839
      bta=parm                                                             1839
      omb=1.0-bta                                                          1839
      dev1=0.0                                                             1839
      dev0=0.0                                                             1840
13360 do 13361 ic=1,nc                                                     1840
      q0=dot_product(w,y(:,ic))                                            1841
      if(q0 .gt. pmin)goto 13381                                           1841
      jerr =8000+ic                                                        1841
      return                                                               1841
13381 continue                                                             1842
      if(q0 .lt. 1.0-pmin)goto 13401                                       1842
      jerr =9000+ic                                                        1842
      return                                                               1842
13401 continue                                                             1843
      if(intr .ne. 0)goto 13421                                            1843
      q0=1.0/nc                                                            1843
      b(0,ic)=0.0                                                          1843
      goto 13431                                                           1844
13421 continue                                                             1844
      b(0,ic)=log(q0)                                                      1844
      dev1=dev1-q0*b(0,ic)                                                 1844
13431 continue                                                             1845
13411 continue                                                             1845
      b(1:ni,ic)=0.0                                                       1846
13361 continue                                                             1847
13362 continue                                                             1847
      if(intr.eq.0) dev1=log(float(nc))                                    1847
      ixx=0                                                                1847
      al=0.0                                                               1848
      if(nonzero(no*nc,g) .ne. 0)goto 13451                                1849
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         1849
      sxp=0.0                                                              1850
13460 do 13461 ic=1,nc                                                     1850
      q(:,ic)=exp(b(0,ic))                                                 1850
      sxp=sxp+q(:,ic)                                                      1850
13461 continue                                                             1851
13462 continue                                                             1851
      goto 13471                                                           1852
13451 continue                                                             1852
13480 do 13481 i=1,no                                                      1852
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1852
13481 continue                                                             1852
13482 continue                                                             1852
      sxp=0.0                                                              1853
      if(intr .ne. 0)goto 13501                                            1853
      b(0,:)=0.0                                                           1853
      goto 13511                                                           1854
13501 continue                                                             1854
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 1854
      if(jerr.ne.0) return                                                 1854
13511 continue                                                             1855
13491 continue                                                             1855
      dev1=0.0                                                             1856
13520 do 13521 ic=1,nc                                                     1856
      q(:,ic)=b(0,ic)+g(:,ic)                                              1857
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             1858
      q(:,ic)=exp(q(:,ic))                                                 1858
      sxp=sxp+q(:,ic)                                                      1859
13521 continue                                                             1860
13522 continue                                                             1860
      sxpl=w*log(sxp)                                                      1860
13530 do 13531 ic=1,nc                                                     1860
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  1860
13531 continue                                                             1861
13532 continue                                                             1861
13471 continue                                                             1862
13441 continue                                                             1862
13540 do 13541 ic=1,nc                                                     1862
13550 do 13551 i=1,no                                                      1862
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               1862
13551 continue                                                             1862
13552 continue                                                             1862
13541 continue                                                             1863
13542 continue                                                             1863
      dev0=dev0+dev1                                                       1864
      if(kopt .le. 0)goto 13571                                            1865
      if(isd .le. 0 .or. intr .eq. 0)goto 13591                            1865
      xv=0.25                                                              1865
      goto 13601                                                           1866
13591 continue                                                             1866
13610 do 13611 j=1,ni                                                      1866
      if(ju(j).ne.0) xv(j,:)=0.25*dot_product(w,x(:,j)**2)                 1866
13611 continue                                                             1866
13612 continue                                                             1866
13601 continue                                                             1867
13581 continue                                                             1867
13571 continue                                                             1869
      alf=1.0                                                              1871
      if(flmin .ge. 1.0)goto 13631                                         1871
      eqs=max(eps,flmin)                                                   1871
      alf=eqs**(1.0/(nlam-1))                                              1871
13631 continue                                                             1872
      m=0                                                                  1872
      mm=0                                                                 1872
      nin=0                                                                1872
      nlp=0                                                                1872
      mnl=min(mnlam,nlam)                                                  1872
      bs=0.0                                                               1872
      shr=shri*dev0                                                        1873
      ga=0.0                                                               1874
13640 do 13641 ic=1,nc                                                     1874
      r=w*(y(:,ic)-q(:,ic)/sxp)                                            1875
13650 do 13651 j=1,ni                                                      1875
      if(ju(j).ne.0) ga(j)=max(ga(j),abs(dot_product(r,x(:,j))))           1875
13651 continue                                                             1876
13652 continue                                                             1876
13641 continue                                                             1877
13642 continue                                                             1877
13660 do 13661 ilm=1,nlam                                                  1877
      al0=al                                                               1878
      if(flmin .lt. 1.0)goto 13681                                         1878
      al=ulam(ilm)                                                         1878
      goto 13671                                                           1879
13681 if(ilm .le. 2)goto 13691                                             1879
      al=al*alf                                                            1879
      goto 13671                                                           1880
13691 if(ilm .ne. 1)goto 13701                                             1880
      al=big                                                               1880
      goto 13711                                                           1881
13701 continue                                                             1881
      al0=0.0                                                              1882
13720 do 13721 j=1,ni                                                      1882
      if(ju(j).eq.0)goto 13721                                             1882
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            1882
13721 continue                                                             1883
13722 continue                                                             1883
      al0=al0/max(bta,1.0d-3)                                              1883
      al=alf*al0                                                           1884
13711 continue                                                             1885
13671 continue                                                             1885
      al2=al*omb                                                           1885
      al1=al*bta                                                           1885
      tlam=bta*(2.0*al-al0)                                                1886
13730 do 13731 k=1,ni                                                      1886
      if(ixx(k).eq.1)goto 13731                                            1886
      if(ju(k).eq.0)goto 13731                                             1887
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     1888
13731 continue                                                             1889
13732 continue                                                             1889
10880 continue                                                             1890
13740 continue                                                             1890
13741 continue                                                             1890
      ix=0                                                                 1890
      jx=ix                                                                1890
      ig=0                                                                 1891
13750 do 13751 ic=1,nc                                                     1891
      bs(0,ic)=b(0,ic)                                                     1892
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          1893
      xmz=0.0                                                              1894
13760 do 13761 i=1,no                                                      1894
      pic=q(i,ic)/sxp(i)                                                   1895
      if(pic .ge. pfm)goto 13781                                           1895
      pic=0.0                                                              1895
      v(i)=0.0                                                             1895
      goto 13771                                                           1896
13781 if(pic .le. pfx)goto 13791                                           1896
      pic=1.0                                                              1896
      v(i)=0.0                                                             1896
      goto 13801                                                           1897
13791 continue                                                             1897
      v(i)=w(i)*pic*(1.0-pic)                                              1897
      xmz=xmz+v(i)                                                         1897
13801 continue                                                             1898
13771 continue                                                             1898
      r(i)=w(i)*(y(i,ic)-pic)                                              1899
13761 continue                                                             1900
13762 continue                                                             1900
      if(xmz.le.vmin)goto 13751                                            1900
      ig=1                                                                 1901
      if(kopt .ne. 0)goto 13821                                            1902
13830 do 13831 j=1,ni                                                      1902
      if(ixx(j).gt.0) xv(j,ic)=dot_product(v,x(:,j)**2)                    1902
13831 continue                                                             1903
13832 continue                                                             1903
13821 continue                                                             1904
13840 continue                                                             1904
13841 continue                                                             1904
      nlp=nlp+1                                                            1904
      dlx=0.0                                                              1905
13850 do 13851 k=1,ni                                                      1905
      if(ixx(k).eq.0)goto 13851                                            1906
      bk=b(k,ic)                                                           1906
      gk=dot_product(r,x(:,k))                                             1907
      u=gk+xv(k,ic)*b(k,ic)                                                1907
      au=abs(u)-vp(k)*al1                                                  1908
      if(au .gt. 0.0)goto 13871                                            1908
      b(k,ic)=0.0                                                          1908
      goto 13881                                                           1909
13871 continue                                                             1910
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))   1912 
     *)
13881 continue                                                             1913
13861 continue                                                             1913
      d=b(k,ic)-bk                                                         1913
      if(abs(d).le.0.0)goto 13851                                          1914
      dlx=max(dlx,xv(k,ic)*d**2)                                           1914
      r=r-d*v*x(:,k)                                                       1915
      if(mm(k) .ne. 0)goto 13901                                           1915
      nin=nin+1                                                            1916
      if(nin .le. nx)goto 13921                                            1916
      jx=1                                                                 1916
      goto 13852                                                           1916
13921 continue                                                             1917
      mm(k)=nin                                                            1917
      m(nin)=k                                                             1918
13901 continue                                                             1919
13851 continue                                                             1920
13852 continue                                                             1920
      if(jx.gt.0)goto 13842                                                1921
      d=0.0                                                                1921
      if(intr.ne.0) d=sum(r)/xmz                                           1922
      if(d .eq. 0.0)goto 13941                                             1922
      b(0,ic)=b(0,ic)+d                                                    1922
      dlx=max(dlx,xmz*d**2)                                                1922
      r=r-d*v                                                              1922
13941 continue                                                             1923
      if(dlx.lt.shr)goto 13842                                             1924
      if(nlp .le. maxit)goto 13961                                         1924
      jerr=-ilm                                                            1924
      return                                                               1924
13961 continue                                                             1925
13970 continue                                                             1925
13971 continue                                                             1925
      nlp=nlp+1                                                            1925
      dlx=0.0                                                              1926
13980 do 13981 l=1,nin                                                     1926
      k=m(l)                                                               1926
      bk=b(k,ic)                                                           1927
      gk=dot_product(r,x(:,k))                                             1928
      u=gk+xv(k,ic)*b(k,ic)                                                1928
      au=abs(u)-vp(k)*al1                                                  1929
      if(au .gt. 0.0)goto 14001                                            1929
      b(k,ic)=0.0                                                          1929
      goto 14011                                                           1930
14001 continue                                                             1931
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))   1933 
     *)
14011 continue                                                             1934
13991 continue                                                             1934
      d=b(k,ic)-bk                                                         1934
      if(abs(d).le.0.0)goto 13981                                          1935
      dlx=max(dlx,xv(k,ic)*d**2)                                           1935
      r=r-d*v*x(:,k)                                                       1936
13981 continue                                                             1937
13982 continue                                                             1937
      d=0.0                                                                1937
      if(intr.ne.0) d=sum(r)/xmz                                           1938
      if(d .eq. 0.0)goto 14031                                             1938
      b(0,ic)=b(0,ic)+d                                                    1939
      dlx=max(dlx,xmz*d**2)                                                1939
      r=r-d*v                                                              1940
14031 continue                                                             1941
      if(dlx.lt.shr)goto 13972                                             1941
      if(nlp .le. maxit)goto 14051                                         1941
      jerr=-ilm                                                            1941
      return                                                               1941
14051 continue                                                             1942
      goto 13971                                                           1943
13972 continue                                                             1943
      goto 13841                                                           1944
13842 continue                                                             1944
      if(jx.gt.0)goto 13752                                                1945
      if(xmz*(b(0,ic)-bs(0,ic))**2.gt.shr) ix=1                            1946
      if(ix .ne. 0)goto 14071                                              1947
14080 do 14081 j=1,nin                                                     1947
      k=m(j)                                                               1948
      if(xv(k,ic)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 14101                1948
      ix=1                                                                 1948
      goto 14082                                                           1948
14101 continue                                                             1949
14081 continue                                                             1950
14082 continue                                                             1950
14071 continue                                                             1951
14110 do 14111 i=1,no                                                      1951
      fi=b(0,ic)+g(i,ic)                                                   1953
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))         1954
      fi=min(max(exmn,fi),exmx)                                            1954
      sxp(i)=sxp(i)-q(i,ic)                                                1955
      q(i,ic)=min(max(emin*sxp(i),exp(fi)),emax*sxp(i))                    1956
      sxp(i)=sxp(i)+q(i,ic)                                                1957
14111 continue                                                             1958
14112 continue                                                             1958
13751 continue                                                             1959
13752 continue                                                             1959
      s=-sum(b(0,:))/nc                                                    1959
      b(0,:)=b(0,:)+s                                                      1959
      di=s                                                                 1960
14120 do 14121 j=1,nin                                                     1960
      l=m(j)                                                               1961
      if(vp(l) .gt. 0.0)goto 14141                                         1961
      s=sum(b(l,:))/nc                                                     1961
      goto 14151                                                           1962
14141 continue                                                             1962
      s=elc(parm,nc,cl(:,l),b(l,:),is)                                     1962
14151 continue                                                             1963
14131 continue                                                             1963
      b(l,:)=b(l,:)-s                                                      1963
      di=di-s*x(:,l)                                                       1964
14121 continue                                                             1965
14122 continue                                                             1965
      di=exp(di)                                                           1965
      sxp=sxp*di                                                           1965
14160 do 14161 ic=1,nc                                                     1965
      q(:,ic)=q(:,ic)*di                                                   1965
14161 continue                                                             1966
14162 continue                                                             1966
      if(jx.gt.0)goto 13742                                                1966
      if(ig.eq.0)goto 13742                                                1967
      if(ix .ne. 0)goto 14181                                              1968
14190 do 14191 k=1,ni                                                      1968
      if(ixx(k).eq.1)goto 14191                                            1968
      if(ju(k).eq.0)goto 14191                                             1968
      ga(k)=0.0                                                            1968
14191 continue                                                             1969
14192 continue                                                             1969
14200 do 14201 ic=1,nc                                                     1969
      r=w*(y(:,ic)-q(:,ic)/sxp)                                            1970
14210 do 14211 k=1,ni                                                      1970
      if(ixx(k).eq.1)goto 14211                                            1970
      if(ju(k).eq.0)goto 14211                                             1971
      ga(k)=max(ga(k),abs(dot_product(r,x(:,k))))                          1972
14211 continue                                                             1973
14212 continue                                                             1973
14201 continue                                                             1974
14202 continue                                                             1974
14220 do 14221 k=1,ni                                                      1974
      if(ixx(k).eq.1)goto 14221                                            1974
      if(ju(k).eq.0)goto 14221                                             1975
      if(ga(k) .le. al1*vp(k))goto 14241                                   1975
      ixx(k)=1                                                             1975
      ix=1                                                                 1975
14241 continue                                                             1976
14221 continue                                                             1977
14222 continue                                                             1977
      if(ix.eq.1) go to 10880                                              1978
      goto 13742                                                           1979
14181 continue                                                             1980
      goto 13741                                                           1981
13742 continue                                                             1981
      if(jx .le. 0)goto 14261                                              1981
      jerr=-10000-ilm                                                      1981
      goto 13662                                                           1981
14261 continue                                                             1981
      devi=0.0                                                             1982
14270 do 14271 ic=1,nc                                                     1983
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          1983
      a0(ic,ilm)=b(0,ic)                                                   1984
14280 do 14281 i=1,no                                                      1984
      if(y(i,ic).le.0.0)goto 14281                                         1985
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           1986
14281 continue                                                             1987
14282 continue                                                             1987
14271 continue                                                             1988
14272 continue                                                             1988
      kin(ilm)=nin                                                         1988
      alm(ilm)=al                                                          1988
      lmu=ilm                                                              1989
      dev(ilm)=(dev1-devi)/dev0                                            1989
      if(ig.eq.0)goto 13662                                                1990
      if(ilm.lt.mnl)goto 13661                                             1990
      if(flmin.ge.1.0)goto 13661                                           1991
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 13662             1992
      if(dev(ilm).gt.devmax)goto 13662                                     1992
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 13662                             1993
13661 continue                                                             1994
13662 continue                                                             1994
      g=log(q)                                                             1994
14290 do 14291 i=1,no                                                      1994
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         1994
14291 continue                                                             1995
14292 continue                                                             1995
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,ga,ixx)                           1996
      return                                                               1997
      end                                                                  1998

      subroutine get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)          772
      implicit double precision(a-h,o-z)                                    773
      data sml0,eps0,big0,mnlam0,rsqmax0,pmin0,exmx0  /1.0d-5,1.0d-6,9.9    775 
     *d35,5,0.999,1.0d-9,250.0/
      sml=sml0                                                              775
      eps=eps0                                                              775
      big=big0                                                              775
      mnlam=mnlam0                                                          775
      rsqmax=rsqmax0                                                        776
      pmin=pmin0                                                            776
      exmx=exmx0                                                            777
      return                                                                778
      entry chg_fract_dev(arg)                                              778
      sml0=arg                                                              778
      return                                                                779
      entry chg_dev_max(arg)                                                779
      rsqmax0=arg                                                           779
      return                                                                780
      entry chg_min_flmin(arg)                                              780
      eps0=arg                                                              780
      return                                                                781
      entry chg_big(arg)                                                    781
      big0=arg                                                              781
      return                                                                782
      entry chg_min_lambdas(irg)                                            782
      mnlam0=irg                                                            782
      return                                                                783
      entry chg_min_null_prob(arg)                                          783
      pmin0=arg                                                             783
      return                                                                784
      entry chg_max_exp(arg)                                                784
      exmx0=arg                                                             784
      return                                                                785
      end                                                                   786
      function dev2(n,w,y,p,pmin)                                          1765
      implicit double precision(a-h,o-z)                                   1766
      double precision w(n),y(n),p(n)                                      1767
      pmax=1.0-pmin                                                        1767
      s=0.0                                                                1768
13340 do 13341 i=1,n                                                       1768
      pi=min(max(pmin,p(i)),pmax)                                          1769
      s=s-w(i)*(y(i)*log(pi)+(1.0-y(i))*log(1.0-pi))                       1770
13341 continue                                                             1771
13342 continue                                                             1771
      dev2=s                                                               1772
      return                                                               1773
      end                                                                  1774
      function azero(n,y,g,q,jerr)                                         1775
      implicit double precision(a-h,o-z)                                   1776
      parameter(eps=1.0d-7)                                                1777
      double precision y(n),g(n),q(n)                                      1778
      double precision, dimension (:), allocatable :: e,p,w                     
      azero = 0.0                                                          1782
      allocate(e(1:n),stat=jerr)                                           1783
      if(jerr.ne.0) return                                                 1784
      allocate(p(1:n),stat=jerr)                                           1785
      if(jerr.ne.0) return                                                 1786
      allocate(w(1:n),stat=jerr)                                           1787
      if(jerr.ne.0) return                                                 1788
      az=0.0                                                               1788
      e=exp(-g)                                                            1788
      qy=dot_product(q,y)                                                  1788
      p=1.0/(1.0+e)                                                        1789
13350 continue                                                             1789
13351 continue                                                             1789
      w=q*p*(1.0-p)                                                        1790
      d=(qy-dot_product(q,p))/sum(w)                                       1790
      az=az+d                                                              1790
      if(abs(d).lt.eps)goto 13352                                          1791
      ea0=exp(-az)                                                         1791
      p=1.0/(1.0+ea0*e)                                                    1792
      goto 13351                                                           1793
13352 continue                                                             1793
      azero=az                                                             1794
      deallocate(e,p,w)                                                    1795
      return                                                               1796
      end                                                                  1797
      function nonzero(n,v)                                                3087
      implicit double precision(a-h,o-z)                                   3088
      double precision v(n)                                                3089
      nonzero=0                                                            3089
18400 do 18401 i=1,n                                                       3089
      if(v(i) .eq. 0.0)goto 18421                                          3089
      nonzero=1                                                            3089
      return                                                               3089
18421 continue                                                             3089
18401 continue                                                             3090
18402 continue                                                             3090
      return                                                               3091
      end                                                                  3092
      subroutine kazero(kk,n,y,g,q,az,jerr)                                1999
      implicit double precision(a-h,o-z)                                   2000
      parameter(eps=1.0d-7)                                                2001
      double precision y(n,kk),g(n,kk),q(n),az(kk)                         2002
      double precision, dimension (:), allocatable :: s                         
      double precision, dimension (:,:), allocatable :: e                       
      allocate(e(1:n,1:kk),stat=jerr)                                           
      if(jerr.ne.0) return                                                      
      allocate(s(1:n),stat=jerr)                                           2009
      if(jerr.ne.0) return                                                 2010
      az=0.0                                                               2010
      e=exp(g)                                                             2010
14300 do 14301 i=1,n                                                       2010
      s(i)=sum(e(i,:))                                                     2010
14301 continue                                                             2011
14302 continue                                                             2011
14310 continue                                                             2011
14311 continue                                                             2011
      dm=0.0                                                               2012
14320 do 14321 k=1,kk                                                      2012
      t=0.0                                                                2012
      u=t                                                                  2013
14330 do 14331 i=1,n                                                       2013
      pik=e(i,k)/s(i)                                                      2014
      t=t+q(i)*(y(i,k)-pik)                                                2014
      u=u+q(i)*pik*(1.0-pik)                                               2015
14331 continue                                                             2016
14332 continue                                                             2016
      d=t/u                                                                2016
      az(k)=az(k)+d                                                        2016
      ed=exp(d)                                                            2016
      dm=max(dm,abs(d))                                                    2017
14340 do 14341 i=1,n                                                       2017
      z=e(i,k)                                                             2017
      e(i,k)=z*ed                                                          2017
      s(i)=s(i)-z+e(i,k)                                                   2017
14341 continue                                                             2018
14342 continue                                                             2018
14321 continue                                                             2019
14322 continue                                                             2019
      if(dm.lt.eps)goto 14312                                              2019
      goto 14311                                                           2020
14312 continue                                                             2020
      az=az-sum(az)/kk                                                     2021
      deallocate(e,s)                                                      2022
      return                                                               2023
      end                                                                  2024
      function elc(parm,n,cl,a,m)                                          2025
      implicit double precision(a-h,o-z)                                   2026
      double precision a(n),cl(2)                                          2026
      integer m(n)                                                         2027
      fn=n                                                                 2027
      am=sum(a)/fn                                                         2028
      if((parm .ne. 0.0) .and. (n .ne. 2))goto 14361                       2028
      elc=am                                                               2028
      go to 14370                                                          2028
14361 continue                                                             2029
14380 do 14381 i=1,n                                                       2029
      m(i)=i                                                               2029
14381 continue                                                             2029
14382 continue                                                             2029
      call psort7(a,m,1,n)                                                 2030
      if(a(m(1)) .ne. a(m(n)))goto 14401                                   2030
      elc=a(1)                                                             2030
      go to 14370                                                          2030
14401 continue                                                             2031
      if(mod(n,2) .ne. 1)goto 14421                                        2031
      ad=a(m(n/2+1))                                                       2031
      goto 14431                                                           2032
14421 continue                                                             2032
      ad=0.5*(a(m(n/2+1))+a(m(n/2)))                                       2032
14431 continue                                                             2033
14411 continue                                                             2033
      if(parm .ne. 1.0)goto 14451                                          2033
      elc=ad                                                               2033
      go to 14370                                                          2033
14451 continue                                                             2034
      b1=min(am,ad)                                                        2034
      b2=max(am,ad)                                                        2034
      k2=1                                                                 2035
14460 continue                                                             2035
14461 if(a(m(k2)).gt.b1)goto 14462                                         2035
      k2=k2+1                                                              2035
      goto 14461                                                           2035
14462 continue                                                             2035
      k1=k2-1                                                              2036
14470 continue                                                             2036
14471 if(a(m(k2)).ge.b2)goto 14472                                         2036
      k2=k2+1                                                              2036
      goto 14471                                                           2037
14472 continue                                                             2037
      r=parm/((1.0-parm)*fn)                                               2037
      is=0                                                                 2037
      sm=n-2*(k1-1)                                                        2038
14480 do 14481 k=k1,k2-1                                                   2038
      sm=sm-2.0                                                            2038
      s=r*sm+am                                                            2039
      if(s .le. a(m(k)) .or. s .gt. a(m(k+1)))goto 14501                   2039
      is=k                                                                 2039
      goto 14482                                                           2039
14501 continue                                                             2040
14481 continue                                                             2041
14482 continue                                                             2041
      if(is .eq. 0)goto 14521                                              2041
      elc=s                                                                2041
      go to 14370                                                          2041
14521 continue                                                             2041
      r2=2.0*r                                                             2041
      s1=a(m(k1))                                                          2041
      am2=2.0*am                                                           2042
      cri=r2*sum(abs(a-s1))+s1*(s1-am2)                                    2042
      elc=s1                                                               2043
14530 do 14531 k=k1+1,k2                                                   2043
      s=a(m(k))                                                            2043
      if(s.eq.s1)goto 14531                                                2044
      c=r2*sum(abs(a-s))+s*(s-am2)                                         2045
      if(c .ge. cri)goto 14551                                             2045
      cri=c                                                                2045
      elc=s                                                                2045
14551 continue                                                             2045
      s1=s                                                                 2046
14531 continue                                                             2047
14532 continue                                                             2047
14370 continue                                                             2047
      elc=max(maxval(a-cl(2)),min(minval(a-cl(1)),elc))                    2048
      return                                                               2049
      end                                                                  2050
      function nintot(ni,nx,nc,a,m,nin,is)                                 2051
      implicit double precision(a-h,o-z)                                   2052
      double precision a(nx,nc)                                            2052
      integer m(nx),is(ni)                                                 2053
      is=0                                                                 2053
      nintot=0                                                             2054
14560 do 14561 ic=1,nc                                                     2054
14570 do 14571 j=1,nin                                                     2054
      k=m(j)                                                               2054
      if(is(k).ne.0)goto 14571                                             2055
      if(a(j,ic).eq.0.0)goto 14571                                         2055
      is(k)=k                                                              2055
      nintot=nintot+1                                                      2056
14571 continue                                                             2056
14572 continue                                                             2056
14561 continue                                                             2057
14562 continue                                                             2057
      return                                                               2058
      end                                                                  2059
      
      subroutine psort7 (v,a,ii,jj)                                             
      implicit double precision(a-h,o-z)                                        
c                                                                               
c     puts into a the permutation vector which sorts v into                     
c     increasing order. the array v is not modified.                            
c     only elements from ii to jj are considered.                               
c     arrays iu(k) and il(k) permit sorting up to 2**(k+1)-1 elements           
c                                                                               
c     this is a modification of cacm algorithm #347 by r. c. singleton,         
c     which is a modified hoare quicksort.                                      
c                                                                               
      dimension a(jj),v(jj),iu(20),il(20)                                       
      integer t,tt                                                              
      integer a                                                                 
      double precision v                                                        
      m=1                                                                       
      i=ii                                                                      
      j=jj                                                                      
 10   if (i.ge.j) go to 80                                                      
 20   k=i                                                                       
      ij=(j+i)/2                                                                
      t=a(ij)                                                                   
      vt=v(t)                                                                   
      if (v(a(i)).le.vt) go to 30                                               
      a(ij)=a(i)                                                                
      a(i)=t                                                                    
      t=a(ij)                                                                   
      vt=v(t)                                                                   
 30   l=j                                                                       
      if (v(a(j)).ge.vt) go to 50                                               
      a(ij)=a(j)                                                                
      a(j)=t                                                                    
      t=a(ij)                                                                   
      vt=v(t)                                                                   
      if (v(a(i)).le.vt) go to 50                                               
      a(ij)=a(i)                                                                
      a(i)=t                                                                    
      t=a(ij)                                                                   
      vt=v(t)                                                                   
      go to 50                                                                  
 40   a(l)=a(k)                                                                 
      a(k)=tt                                                                   
 50   l=l-1                                                                     
      if (v(a(l)).gt.vt) go to 50                                               
      tt=a(l)                                                                   
      vtt=v(tt)                                                                 
 60   k=k+1                                                                     
      if (v(a(k)).lt.vt) go to 60                                               
      if (k.le.l) go to 40                                                      
      if (l-i.le.j-k) go to 70                                                  
      il(m)=i                                                                   
      iu(m)=l                                                                   
      i=k                                                                       
      m=m+1                                                                     
      go to 90                                                                  
 70   il(m)=k                                                                   
      iu(m)=j                                                                   
      j=l                                                                       
      m=m+1                                                                     
      go to 90                                                                  
 80   m=m-1                                                                     
      if (m.eq.0) return                                                        
      i=il(m)                                                                   
      j=iu(m)                                                                   
 90   if (j-i.gt.10) go to 20                                                   
      if (i.eq.ii) go to 10                                                     
      i=i-1                                                                     
 100  i=i+1                                                                     
      if (i.eq.j) go to 80                                                      
      t=a(i+1)                                                                  
      vt=v(t)                                                                   
      if (v(a(i)).le.vt) go to 100                                              
      k=i                                                                       
 110  a(k+1)=a(k)                                                               
      k=k-1                                                                     
      if (vt.lt.v(a(k))) go to 110                                              
      a(k+1)=t                                                                  
      go to 100                                                                 
      end       

! fishnet starts here ---------------------------------------------------------------
      subroutine pfishnet (parm,no,ni,x,g,w,jd,vp,cl,ne,nx,nlam,flmin,u    2918 
     *lam,thr,  isd, maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,ncp,pp,
     *jerr)
      implicit double precision(a-h,o-z)                                   2919
      double precision x(no,ni),y(ncp),g(ncp),w(ncp),vp(ni),ulam(nlam)     2920
      double precision ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam),cl(2,ni)   2921
      integer jd(*),ia(nx),nin(nlam)                                       2922
      double precision, dimension (:), allocatable :: xm,xs,ww,vq               
      integer, dimension (:), allocatable :: ju
      
      integer intr
      double precision xh(ncp,ni)
      integer pp(ncp,2)
      intr = 0 
      y = 0.0
      ! print *, "flmin = ", flmin
      xh =  x(pp(:,1),:) - x(pp(:,2),:)

      if(maxval(vp) .gt. 0.0)goto 17701                                    2926
      jerr=10000                                                           2926
      return                                                               2926
17701 continue                                                             2927
      if(minval(y) .ge. 0.0)goto 17721                                     2927
      jerr=8888                                                            2927
      return                                                               2927
17721 continue                                                             2928
      allocate(ww(1:ncp),stat=jerr)                                        2929
      if(jerr.ne.0) return                                                 2930
      allocate(ju(1:ni),stat=jerr)                                         2931
      if(jerr.ne.0) return                                                 2932
      allocate(vq(1:ni),stat=jerr)                                         2933
      if(jerr.ne.0) return                                                 2934
      allocate(xm(1:ni),stat=jerr)                                         2935
      if(jerr.ne.0) return                                                 2936
      if(isd .le. 0)goto 17741                                             2936
      allocate(xs(1:ni),stat=jerr)                                         2936
      if(jerr.ne.0) return                                                 2936
17741 continue                                                             2937
      call chkvars(ncp,ni,xh,ju)                                           2938
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 2939
      if(maxval(ju) .gt. 0)goto 17761                                      2939
      jerr=7777                                                            2939
      go to 12180                                                          2939
17761 continue                                                             2940
      vq=max(0d0,vp)                                                       2940
      vq=vq*ni/sum(vq)                                                     2941
      ww=max(0d0,w)                                                        2941
      sw=sum(ww)                                                           2941
      if(sw .gt. 0.0)goto 17781                                            2941
      jerr=9999                                                            2941
      go to 12180                                                          2941
17781 continue                                                             2942
      ww=ww/sw                                                             2943
      call lstandard1(ncp,ni,xh,ww,ju,isd,intr,xm,xs)                      2944
      if(isd .le. 0)goto 17801                                             2944
17810 do 17811 j=1,ni                                                      2944
      cl(:,j)=cl(:,j)*xs(j)                                                2944
17811 continue                                                             2944
17812 continue                                                             2944
17801 continue                                                             2945
      call fishnet1(parm,ncp,ni,xh,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam   2947 
     *,thr,  isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) go to 12180                                            2947
      dev0=2.0*sw*dev0                                                     2948
17820 do 17821 k=1,lmu                                                     2948
      nk=nin(k)                                                            2949
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      2950
      if(intr .ne. 0)goto 17841                                            2950
      a0(k)=0.0                                                            2950
      goto 17851                                                           2951
17841 continue                                                             2951
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     2951
17851 continue                                                             2952
17831 continue                                                             2952
17821 continue                                                             2953
17822 continue                                                             2953
12180 continue                                                             2953
      deallocate(ww,ju,vq,xm)                                              2953
      if(isd.gt.0) deallocate(xs)                                          2954
      return                                                               2955
      end                                                                  2956

      subroutine fishnet1(parm,no,ni,x,y,g,q,ju,vp,cl,ne,nx,nlam,flmin,u   2958 
     *lam,shri,  isd,intr,maxit,lmu,a0,ca,m,kin,dev0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   2959
      double precision x(no,ni),y(no),g(no),q(no),vp(ni),ulam(nlam)        2960
      double precision ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam),cl(2,ni)   2961
      integer ju(ni),m(nx),kin(nlam)                                       2962
      double precision, dimension (:), allocatable :: t,w,wr,v,a,f,as,ga        
      integer, dimension (:), allocatable :: mm,ixx                             
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               2966
      sml=sml*10.0                                                         2967
      allocate(a(1:ni),stat=jerr)                                          2968
      if(jerr.ne.0) return                                                 2969
      allocate(as(1:ni),stat=jerr)                                         2970
      if(jerr.ne.0) return                                                 2971
      allocate(t(1:no),stat=jerr)                                          2972
      if(jerr.ne.0) return                                                 2973
      allocate(mm(1:ni),stat=jerr)                                         2974
      if(jerr.ne.0) return                                                 2975
      allocate(ga(1:ni),stat=jerr)                                         2976
      if(jerr.ne.0) return                                                 2977
      allocate(ixx(1:ni),stat=jerr)                                        2978
      if(jerr.ne.0) return                                                 2979
      allocate(wr(1:no),stat=jerr)                                         2980
      if(jerr.ne.0) return                                                 2981
      allocate(v(1:ni),stat=jerr)                                          2982
      if(jerr.ne.0) return                                                 2983
      allocate(w(1:no),stat=jerr)                                          2984
      if(jerr.ne.0) return                                                 2985
      allocate(f(1:no),stat=jerr)                                          2986
      if(jerr.ne.0) return                                                 2987
      bta=parm                                                             2987
      omb=1.0-bta                                                          2988
      t=q*y                                                                2988
      yb=sum(t)                                                            2988
      fmax=log(huge(bta)*0.1)                                              2989
      if(nonzero(no,g) .ne. 0)goto 17871                                   2990
      if(intr .eq. 0)goto 17891                                            2990
      w=q*yb                                                               2990
      az=log(yb)                                                           2990
      f=az                                                                 2990
      dv0=yb*(az-1.0)                                                      2990
      goto 17901                                                           2991
17891 continue                                                             2991
      w=q                                                                  2991
      az=0.0                                                               2991
      f=az                                                                 2991
      dv0=-1.0                                                             2991
17901 continue                                                             2992
17881 continue                                                             2992
      goto 17911                                                           2993
17871 continue                                                             2993
      w=q*exp(sign(min(abs(g),fmax),g))                                    2993
      v0=sum(w)                                                            2994
      if(intr .eq. 0)goto 17931                                            2994
      eaz=yb/v0                                                            2994
      w=eaz*w                                                              2994
      az=log(eaz)                                                          2994
      f=az+g                                                               2995
      dv0=dot_product(t,g)-yb*(1.0-az)                                     2996
      goto 17941                                                           2997
17931 continue                                                             2997
      az=0.0                                                               2997
      f=g                                                                  2997
      dv0=dot_product(t,g)-v0                                              2997
17941 continue                                                             2998
17921 continue                                                             2998
17911 continue                                                             2999
17861 continue                                                             2999
      a=0.0                                                                2999
      as=0.0                                                               2999
      wr=t-w                                                               2999
      v0=1.0                                                               2999
      if(intr.ne.0) v0=yb                                                  2999
      dvr=-yb                                                              3000
17950 do 17951 i=1,no                                                      3000
      if(t(i).gt.0.0) dvr=dvr+t(i)*log(y(i))                               3000
17951 continue                                                             3000
17952 continue                                                             3000
      dvr=dvr-dv0                                                          3000
      dev0=dvr                                                             3002
      alf=1.0                                                              3004
      if(flmin .ge. 1.0)goto 17971                                         3004
      eqs=max(eps,flmin)                                                   3004
      alf=eqs**(1.0/(nlam-1))                                              3004
17971 continue                                                             3005
      m=0                                                                  3005
      mm=0                                                                 3005
      nlp=0                                                                3005
      nin=nlp                                                              3005
      mnl=min(mnlam,nlam)                                                  3005
      shr=shri*dev0                                                        3005
      ixx=0                                                                3005
      al=0.0                                                               3006
17980 do 17981 j=1,ni                                                      3006
      if(ju(j).eq.0)goto 17981                                             3006
      ga(j)=abs(dot_product(wr,x(:,j)))                                    3006
17981 continue                                                             3007
17982 continue                                                             3007
17990 do 17991 ilm=1,nlam                                                  3007
      al0=al                                                               3008
      if(flmin .lt. 1.0)goto 18011                                         3008
      al=ulam(ilm)                                                         3008
      goto 18001                                                           3009
18011 if(ilm .le. 2)goto 18021                                             3009
      al=al*alf                                                            3009
      goto 18001                                                           3010
18021 if(ilm .ne. 1)goto 18031                                             3010
      al=big                                                               3010
      goto 18041                                                           3011
18031 continue                                                             3011
      al0=0.0                                                              3012
18050 do 18051 j=1,ni                                                      3012
      if(ju(j).eq.0)goto 18051                                             3012
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            3012
18051 continue                                                             3013
18052 continue                                                             3013
      al0=al0/max(bta,1.0d-3)                                              3013
      al=alf*al0                                                           3014
18041 continue                                                             3015
18001 continue                                                             3015
      al2=al*omb                                                           3015
      al1=al*bta                                                           3015
      tlam=bta*(2.0*al-al0)                                                3016
18060 do 18061 k=1,ni                                                      3016
      if(ixx(k).eq.1)goto 18061                                            3016
      if(ju(k).eq.0)goto 18061                                             3017
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     3018
18061 continue                                                             3019
18062 continue                                                             3019
10880 continue                                                             3020
18070 continue                                                             3020
18071 continue                                                             3020
      az0=az                                                               3021
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                3022
18080 do 18081 j=1,ni                                                      3022
      if(ixx(j).ne.0) v(j)=dot_product(w,x(:,j)**2)                        3022
18081 continue                                                             3023
18082 continue                                                             3023
18090 continue                                                             3023
18091 continue                                                             3023
      nlp=nlp+1                                                            3023
      dlx=0.0                                                              3024
18100 do 18101 k=1,ni                                                      3024
      if(ixx(k).eq.0)goto 18101                                            3024
      ak=a(k)                                                              3025
      u=dot_product(wr,x(:,k))+v(k)*ak                                     3025
      au=abs(u)-vp(k)*al1                                                  3026
      if(au .gt. 0.0)goto 18121                                            3026
      a(k)=0.0                                                             3026
      goto 18131                                                           3027
18121 continue                                                             3028
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))           3029
18131 continue                                                             3030
18111 continue                                                             3030
      if(a(k).eq.ak)goto 18101                                             3030
      d=a(k)-ak                                                            3030
      dlx=max(dlx,v(k)*d**2)                                               3031
      wr=wr-d*w*x(:,k)                                                     3031
      f=f+d*x(:,k)                                                         3032
      if(mm(k) .ne. 0)goto 18151                                           3032
      nin=nin+1                                                            3032
      if(nin.gt.nx)goto 18102                                              3033
      mm(k)=nin                                                            3033
      m(nin)=k                                                             3034
18151 continue                                                             3035
18101 continue                                                             3036
18102 continue                                                             3036
      if(nin.gt.nx)goto 18092                                              3037
      if(intr .eq. 0)goto 18171                                            3037
      d=sum(wr)/v0                                                         3038
      az=az+d                                                              3038
      dlx=max(dlx,v0*d**2)                                                 3038
      wr=wr-d*w                                                            3038
      f=f+d                                                                3039
18171 continue                                                             3040
      if(dlx.lt.shr)goto 18092                                             3040
      if(nlp .le. maxit)goto 18191                                         3040
      jerr=-ilm                                                            3040
      return                                                               3040
18191 continue                                                             3041
18200 continue                                                             3041
18201 continue                                                             3041
      nlp=nlp+1                                                            3041
      dlx=0.0                                                              3042
18210 do 18211 l=1,nin                                                     3042
      k=m(l)                                                               3042
      ak=a(k)                                                              3043
      u=dot_product(wr,x(:,k))+v(k)*ak                                     3043
      au=abs(u)-vp(k)*al1                                                  3044
      if(au .gt. 0.0)goto 18231                                            3044
      a(k)=0.0                                                             3044
      goto 18241                                                           3045
18231 continue                                                             3046
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))           3047
18241 continue                                                             3048
18221 continue                                                             3048
      if(a(k).eq.ak)goto 18211                                             3048
      d=a(k)-ak                                                            3048
      dlx=max(dlx,v(k)*d**2)                                               3049
      wr=wr-d*w*x(:,k)                                                     3049
      f=f+d*x(:,k)                                                         3051
18211 continue                                                             3051
18212 continue                                                             3051
      if(intr .eq. 0)goto 18261                                            3051
      d=sum(wr)/v0                                                         3051
      az=az+d                                                              3052
      dlx=max(dlx,v0*d**2)                                                 3052
      wr=wr-d*w                                                            3052
      f=f+d                                                                3053
18261 continue                                                             3054
      if(dlx.lt.shr)goto 18202                                             3054
      if(nlp .le. maxit)goto 18281                                         3054
      jerr=-ilm                                                            3054
      return                                                               3054
18281 continue                                                             3055
      goto 18201                                                           3056
18202 continue                                                             3056
      goto 18091                                                           3057
18092 continue                                                             3057
      if(nin.gt.nx)goto 18072                                              3058
      w=q*exp(sign(min(abs(f),fmax),f))                                    3058
      v0=sum(w)                                                            3058
      wr=t-w                                                               3059
      if(v0*(az-az0)**2 .ge. shr)goto 18301                                3059
      ix=0                                                                 3060
18310 do 18311 j=1,nin                                                     3060
      k=m(j)                                                               3061
      if(v(k)*(a(k)-as(k))**2.lt.shr)goto 18311                            3061
      ix=1                                                                 3061
      goto 18312                                                           3062
18311 continue                                                             3063
18312 continue                                                             3063
      if(ix .ne. 0)goto 18331                                              3064
18340 do 18341 k=1,ni                                                      3064
      if(ixx(k).eq.1)goto 18341                                            3064
      if(ju(k).eq.0)goto 18341                                             3065
      ga(k)=abs(dot_product(wr,x(:,k)))                                    3066
      if(ga(k) .le. al1*vp(k))goto 18361                                   3066
      ixx(k)=1                                                             3066
      ix=1                                                                 3066
18361 continue                                                             3067
18341 continue                                                             3068
18342 continue                                                             3068
      if(ix.eq.1) go to 10880                                              3069
      goto 18072                                                           3070
18331 continue                                                             3071
18301 continue                                                             3072
      goto 18071                                                           3073
18072 continue                                                             3073
      if(nin .le. nx)goto 18381                                            3073
      jerr=-10000-ilm                                                      3073
      goto 17992                                                           3073
18381 continue                                                             3074
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               3074
      kin(ilm)=nin                                                         3075
      a0(ilm)=az                                                           3075
      alm(ilm)=al                                                          3075
      lmu=ilm                                                              3076
      dev(ilm)=(dot_product(t,f)-v0-dv0)/dvr                               3077
      if(ilm.lt.mnl)goto 17991                                             3077
      if(flmin.ge.1.0)goto 17991                                           3078
      me=0                                                                 3078
18390 do 18391 j=1,nin                                                     3078
      if(ca(j,ilm).ne.0.0) me=me+1                                         3078
18391 continue                                                             3078
18392 continue                                                             3078
      if(me.gt.ne)goto 17992                                               3079
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 17992              3080
      if(dev(ilm).gt.devmax)goto 17992                                     3081
17991 continue                                                             3082
17992 continue                                                             3082
      g=f                                                                  3083
12180 continue                                                             3083
      deallocate(t,w,wr,v,a,f,as,mm,ga,ixx)                                3084
      return                                                               3085
      end                                                                  3086





      