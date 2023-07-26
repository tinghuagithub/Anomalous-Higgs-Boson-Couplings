
       subroutine ten_red2_forGram(p1sq,B012,Bij)
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c************************************************************************************
       implicit none
       Real*8 p1sq,p1sq2,p1sq3,p1sq4,p1sq5
       Complex*16 B012
       Complex*16 Bij(6,11)
       Real*8 Inv2,Inv5,Inv12,Inv18,Inv36,Inv72,Inv300,Inv600,Inv3600
       Real*8 Inv360,Inv7200,Inv980,Inv35280,Inv352800,Inv705600
       Real*8 Inv560,Inv7840, Inv141120,Inv1411200
       Real*8 Inv22680,Inv90720,Inv1270080,Inv7620480,Inv76204800
       Real*8 Inv25200, Inv453600,Inv1814400,Inv25401600,Inv152409600
       Real*8 Inv304920, Inv6098400, Inv109771200,Inv439084800
       Real*8 Inv6147187200,Inv36883123200,Inv332640,Inv7318080
       real*8 Inv146361600,Inv2634508800,Inv5269017600,Inv73766246400

       If(abs(p1sq).gt.1d-5) then
        Inv2=1d0/2d0
        Inv5=1d0/5d0
        Inv12=1d0/12d0
        Inv18=1d0/18d0
        Inv36=1d0/36d0
        Inv72=1d0/72d0
        Inv300=1d0/300d0
        Inv360=1d0/360d0
        Inv600=1d0/600d0
        Inv3600=1d0/3600d0
        Inv7200=1d0/7200d0
        Inv980=1d0/980d0
        Inv35280=1d0/35280d0
        Inv352800=1d0/352800d0
        Inv705600=1d0/705600d0
        Inv560=1d0/560d0
        Inv7840=1d0/7840d0
        Inv141120=1d0/141120d0
        Inv1411200=1d0/1411200d0
        Inv22680=1d0/22680d0
        Inv90720=1d0/90720d0
        Inv1270080=1d0/1270080d0
        Inv7620480=1d0/7620480d0
        Inv76204800=1d0/76204800d0
        
        Inv25200=1d0/25200d0
        Inv453600=1d0/453600d0
        Inv1814400  =1D0/1814400d0
        Inv25401600  =1D0/25401600d0
        Inv152409600  =1D0/152409600d0
        Inv304920  =1D0/304920d0
        Inv6098400  =1D0/6098400d0
        Inv109771200  =1D0/109771200d0
                         
        Inv439084800    =    1d0/439084800d0
        Inv6147187200    =    1d0/6147187200d0
        Inv36883123200    =    1d0/36883123200d0
        Inv332640    =    1d0/332640d0
        Inv7318080    =    1d0/7318080d0
        Inv146361600    =    1d0/146361600d0
        Inv2634508800    =    1d0/2634508800d0
        Inv5269017600    =    1d0/5269017600d0
        Inv73766246400    =    1d0/73766246400d0


        p1sq2=p1sq*p1sq
        p1sq3=p1sq2*p1sq
        p1sq4=p1sq3*p1sq
        p1sq5=p1sq4*p1sq

       Bij(1,1)=-(B012*Inv2)
       Bij(1,2)=(1d0+6d0*B012)*Inv18
       Bij(2,2)=-((2d0+3d0*B012)*Inv36*p1sq)
       Bij(1,3)=(-1d0-3d0*B012)*Inv12
       Bij(2,3)=(2d0+3d0*B012)*Inv72*p1sq
       Bij(1,4)=29d0*Inv300+B012*Inv5
       Bij(2,4)=-((11d0+15d0*B012)*Inv600*p1sq)
       Bij(3,4)=(16d0+15d0*B012)*Inv3600*p1sq2
       Bij(1,5)=-Inv360*(60d0*B012+37d0)
       Bij(2,5)=Inv3600*(49d0+60d0*B012)*p1sq
       Bij(3,5)=-Inv7200*(16d0+15d0*B012)*p1sq2
       Bij(1,6)= Inv980*(103d0+140d0*B012)
       Bij(2,6)=-Inv35280*(379d0+420d0*B012)*p1sq
       Bij(3,6)=Inv352800*(463d0+420d0*B012)*p1sq2
       Bij(4,6)=-Inv705600*(142d0+105d0*B012)*p1sq3 
       Bij(1,7)=-Inv560*(59d0+70d0*B012)
       Bij(2,7)=Inv7840*(69d0+70d0*B012)*p1sq
       Bij(3,7)=-Inv141120*(121d0+105d0*B012)*p1sq2
       Bij(4,7)=Inv1411200*(142d0+105d0*B012)*p1sq3
       Bij(1,8)=Inv22680*(2369d0+2520d0*B012)
       Bij(2,8)=-Inv90720*(671d0+630d0*B012)*p1sq
       Bij(3,8)=Inv1270080*(761d0+630d0*B012)*p1sq2
       Bij(4,8)=-Inv7620480*(433d0+315d0*B012)*p1sq3
       Bij(5,8)=Inv76204800*(496d0+315d0*B012)*p1sq4

       Bij(1,9)=(-2593d0 - 2520d0*B012)*Inv25200
       Bij(2,9)=(2873d0 + 2520d0*B012)*Inv453600*p1sq
       Bij(3,9)=-((797d0 + 630d0*B012)*Inv1814400*p1sq2)
       Bij(4,9)=(887d0 + 630d0*B012)*Inv25401600*p1sq3
       Bij(5,9)=-((496d0 + 315d0*B012)*Inv152409600*p1sq4)
                 
       Bij(1,10)=(30791d0 + 27720d0*B012)*Inv304920
       Bij(2,10)=-((33563d0 + 27720d0*B012)*Inv6098400*p1sq)
       Bij(3,10)=(36643d0 + 27720d0*B012)*Inv109771200*p1sq2
       Bij(4,10)=-((10027d0 + 6930d0*B012)*Inv439084800*p1sq3)
       Bij(5,10)=(11017d0 + 6930d0*B012)*Inv6147187200*p1sq4
       Bij(6,10)=-((6086d0 + 3465d0*B012)*Inv36883123200*p1sq5)
                 
       Bij(1,11)=(-32891d0 - 27720d0*B012)*Inv332640
       Bij(2,11)=(35411d0 + 27720d0*B012)*Inv7318080*p1sq
       Bij(3,11)=-((38183d0 + 27720d0*B012)*Inv146361600*p1sq2)
       Bij(4,11)=(41263d0 + 27720d0*B012)*Inv2634508800*p1sq3
       Bij(5,11)=-((5591d0 + 3465d0*B012)*Inv5269017600*p1sq4)
       Bij(6,11)=(6086d0 + 3465d0*B012)*Inv73766246400*p1sq5

       return
       else
       Bij(1,1)=0d0
       Bij(1,2)=0d0
       Bij(2,2)=0d0
       Bij(1,3)=0d0
       Bij(2,3)=0d0
       Bij(1,4)=0d0
       Bij(2,4)=0d0
       Bij(3,4)=0d0
       Bij(1,5)=0d0
       Bij(2,5)=0d0
       Bij(3,5)=0d0
       Bij(1,6)=0d0
       Bij(2,6)=0d0
       Bij(3,6)=0d0
       Bij(4,6)=0d0
       Bij(1,7)=0d0
       Bij(2,7)=0d0
       Bij(3,7)=0d0
       Bij(4,7)=0d0
       Bij(1,8)=0d0
       Bij(2,8)=0d0
       Bij(3,8)=0d0
       Bij(4,8)=0d0
       Bij(5,8)=0d0

       Bij(1,9)=0d0
       Bij(2,9)=0d0
       Bij(3,9)=0d0
       Bij(4,9)=0d0
       Bij(5,9)=0d0
                 
       Bij(1,10)=0d0
       Bij(2,10)=0d0
       Bij(3,10)=0d0
       Bij(4,10)=0d0
       Bij(5,10)=0d0
       Bij(6,10)=0d0
                 
       Bij(1,11)=0d0
       Bij(2,11)=0d0
       Bij(3,11)=0d0
       Bij(4,11)=0d0
       Bij(5,11)=0d0
       Bij(6,11)=0d0

       Endif
       return
       End
