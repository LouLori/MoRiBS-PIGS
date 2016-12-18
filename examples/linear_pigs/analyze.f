      implicit real*8(a-h,o-z)
      parameter (nn=nn_input)
      dimension pot(nn)
      dimension energy(nn)
      open(1,file='analyze_input')
      open(2,file='average_analyze_input')
      param1 = param2
      sum1=0.0d0
      sum1var=0.0d0
      sum2=0.0d0
      sum2var=0.0d0
      do i=1,nn
        read(1,*)iblock,ake,pot(i),energy(i),aa,bb,cc,dd,ee,ff
        sum1=sum1+pot(i)
        sum1var=sum1var+pot(i)*pot(i)
        sum2=sum2+energy(i)
        sum2var=sum2var+energy(i)*energy(i)
      enddo
      sum11=sum1/float(nn)
      sum11var=sum1var/float(nn)
      sum22=sum2/float(nn)
      sum22var=sum2var/float(nn)
      dev1=sqrt(sum11var-sum11*sum11)
      dev2=sqrt(sum22var-sum22*sum22)
      write(2,*)nn,param1,sum11,sum22,dev1/float(nn),dev2/float(nn)
      stop
      end
