!This program creates an optimized WS-DD sequence for 10 multichannel instruments such as the Syscal Pro

program main
    implicit none

    !Define variables

    INTEGER::channels,roll_ch,opt_ch,max_ap,ii,i,qp_C
    INTEGER,dimension(100000,2)::ab
    INTEGER,dimension(100000,4)::abmn,abmnO,abmnR
    INTEGER::ioEL,ioElr,roll_count,io

!Ask the input values to the user


print *, 'This program produces an optimized DD+WS sequence for ERT surveys:'

  print *, 'How many channels you have?:'
  read(*,*) channels

  print *, 'How many channels you roll along?:'
  read(*,*) roll_ch

  print *, 'How many channels does your instrument measure at the same time?:'
  read(*,*) opt_ch


    !Define the aperture vector
call apertures(channels,ab,max_ap)


    !Calculate the quadripoles for each aperture
call quadripoles (ab,abmn,max_ap,channels,opt_ch,qp_C)


    !Make it optimized for removing fast the first cable
call reorderEle(abmn,qp_C,abmnO,channels) !abmnO is the output of this subroutine
abmn(1:qp_C,:)=abmnO(1:qp_C,:)


    !Create the roll along sequence
call rollalong(abmn,roll_ch,abmnR,qp_C,roll_count)


    !Save the results in a format for ElectrePro

open(newunit=ioEl, file="DDWSopt.txt", status="new", action="write")

write(ioEl,fmt='(G0.1)') '# X Y Z'

do i =1,channels
write( ioEl,fmt='(*(G0.1,:,Tr1,:))') i, i-1,0,0


end do


write(ioEl,fmt='(G0.1)') '# A B M N'

do i =1,qp_C-1

write( ioEl,fmt='(*(G0.1,:,Tr1,:))') i,abmn(i,1),abmn(i,2),abmn(i,3),abmn(i,4)

end do

close(ioEl)

!Save roll along

open(newunit=ioElr, file="DDWSoptR.txt", status="new", action="write")

write(ioElr,fmt='(G0.1)') '# X Y Z'

do i =1,channels
write( ioElr,fmt='(*(G0.1,:,Tr1,:))') i, i-1,0,0


end do

print*,roll_count
write(ioElr,fmt='(G0.1)') '# A B M N'

do i =1,roll_count-1

write( ioElr,fmt='(*(G0.1,:,Tr1,:))') i,abmnR(i,1),abmnR(i,2),abmnR(i,3),abmnR(i,4)

end do

close(ioElr)



!save csv formatted for resipy forward.


open(newunit=io, file="DDWSopt.csv", status="new", action="write")

write(io,fmt='(G0.1)') 'A,B,M,N'

do i =1,qp_C
write( io,fmt='(*(G0.1,:,",",:))') abmn(i,1),abmn(i,2),abmn(i,3),abmn(i,4)
end do

close(io)




end program main




SUBROUTINE rollalong (abmn,roll_ch,abmnR,qp_C,roll_count)
    integer, dimension(100000,4), intent (in) :: abmn ! input
    integer, dimension(100000,4), intent (out) :: abmnR ! output
    integer,intent(out)::roll_count
    integer,intent(in)::qp_C,roll_ch

INTEGER,dimension(100000)::A,B,M,N,Ar,Br,Mr,Nr
INTEGER::i,j,is_present

meas_count=qp_C
roll_count=0

!compute the roll along sequence, equal to the normal sequence but without those measurements already performed in the previous non-roll measurement.


A=abmn(:,1)
B=abmn(:,2)
M=abmn(:,3)
N=abmn(:,4)

outer:do i=1,meas_count
is_present=0
inner:do j=1,meas_count

if(i/=j .and. A(i)+roll_ch==A(j) .and. B(i)+roll_ch==B(j) .and. M(i)+roll_ch==M(j) .and. N(i)+roll_ch==N(j)) then
is_present=1
end if

end do inner
if (is_present==0) then
roll_count=roll_count+1

Ar(roll_count)=A(i)
Br(roll_count)=B(i)
Mr(roll_count)=M(i)
Nr(roll_count)=N(i)
end if
end do outer

abmnR(:,1)=Ar
abmnR(:,2)=Br
abmnR(:,3)=Mr
abmnR(:,4)=Nr

END SUBROUTINE




SUBROUTINE reorderEle(abmn,qp_C,abmnO,channels)
    integer, dimension(100000,4), intent (in) :: abmn ! input
    integer, dimension(100000,4), intent (out) :: abmnO ! output
    integer,intent(in)::qp_C,channels

INTEGER::groupcount,groupminA,groupminB,groupminM,groupminN,indmin,indmax,i,meas_count
INTEGER,dimension(100000)::group,groupmin,A,B,M,N,Aord,Bord,Mord,Nord


meas_count=qp_C
A=abmn(:,1)
B=abmn(:,2)
M=abmn(:,3)
N=abmn(:,4)



! order the sequence in a way one can remove the first cable the sooner possible

groupcount=1
group(1)=1

! calculate the first index of each optimized group of measurements

do i=2, meas_count

if(A(i)/=A(i-1) .or.B(i)/=B(i-1)) then


groupcount=groupcount+1
group(groupcount)=i
end if
end do

! group 1 quindi va da 1 a group(2)-1, e così via; group finale è group(groupcount):meas_count
!di ciascun gruppo devo calcolare il numero più basso

!Calculate the minimum electrode index for each group of measurements

group(groupcount+1)=meas_count+1

do i=1,groupcount

groupminA=minval(A(group(i):(group(i+1)-1)))
groupminB=minval(B(group(i):(group(i+1)-1)))
groupminM=minval(M(group(i):(group(i+1)-1)))
groupminN=minval(N(group(i):(group(i+1)-1)))
groupmin(i)=min(groupminA,groupminB,groupminM,groupminN)

end do


!re-write the A, B, M, N arrays writing first those that have the least minimum electrode index.
!put the results in Aord, Bord, Mord, Nord arrays

indmin=1
indmax=0

do i=1, channels

do j=1,groupcount
if(groupmin(j)==i) then

indmax=indmin+group(j+1)-group(j)-1

Aord(indmin:indmax)=A(group(j):(group(j+1)-1))
Bord(indmin:indmax)=B(group(j):(group(j+1)-1))
Mord(indmin:indmax)=M(group(j):(group(j+1)-1))
Nord(indmin:indmax)=N(group(j):(group(j+1)-1))

indmin=indmax+1

end if
end do
end do

abmnO(:,1)=Aord
abmnO(:,2)=Bord

abmnO(:,3)=Mord
abmnO(:,4)=Nord


END SUBROUTINE


SUBROUTINE quadripoles(ab,abmn,max_ap,channels,opt_ch,qp_C)
    integer, dimension(100000,2), intent (in) :: ab ! input
    integer,  intent (in) :: max_ap,channels,opt_ch ! input
    integer, dimension(100000,4), intent (out) :: abmn ! output
    integer,intent(out) :: qp_C !output
    INTEGER::qp_I,ab_ap,k,qp_j,kk
    k=floor(channels/(opt_ch-3+1.0)) ! K is a factor above which it is best to perform just WS and not WS+DD. Vedi quadernetto con la tartaruga per calcolo di k

    opt_ch_q=opt_ch

    if(channels<36 .and. opt_ch>5)then  !This will set that, for small number of channels, the algorithm accepts also not fully optimized (10) sequences, to add a little more quadripoles.
     opt_ch_q=floor(channels/10+1.0)!check later where the code uses opt_ch (the maximum) or opt_ch_q, which can be 5,4,3...
      end if

    qp_C=1

    DO qp_I=1,max_ap  !questi sono gli AB in input

    ab_ap=abs(ab(qp_i,2)-ab(qp_i,1))
    if (ab_ap<3) then  !This represents the case in which a small aperture of AB is treated as just DD sequence (WS impossible)

        if(ab(qp_I,2)<=(channels-(opt_ch_q+1)*ab_ap)) then
        do qp_j=1,opt_ch

            abmn(qp_C,1)=ab(qp_I,1)
            abmn(qp_C,2)=ab(qp_I,2)
            abmn(qp_C,3)=ab(qp_I,2)+ab_ap*qp_j
            abmn(qp_C,4)=ab(qp_I,2)+ab_ap*qp_j+ab_ap
            if (abmn(qp_C,4)<=channels) then
            qp_C=qp_C+1
            end if

        end do
        end if

         if(ab(qp_I,1)>((opt_ch_q+1)*ab_ap)) then
        do qp_j=1,opt_ch

            abmn(qp_C,1)=ab(qp_I,1)
            abmn(qp_C,2)=ab(qp_I,2)
            abmn(qp_C,3)=ab(qp_I,1)-ab_ap*qp_j
            abmn(qp_C,4)=ab(qp_I,1)-ab_ap*qp_j-ab_ap
           if (abmn(qp_C,4)>0) then
           qp_C=qp_C+1
           end if

        end do
        end if


   end if

   if (ab_ap>=3 .and. ab_ap<=k) then


   if(ab(qp_I,2)<=(channels-(opt_ch_q-2+1)*ab_ap)) then

        abmn(qp_C,1)=ab(qp_I,1)
        abmn(qp_C,2)=ab(qp_I,2)
        abmn(qp_C,3)=ab(qp_I,1)+1
        abmn(qp_C,4)=ab(qp_I,2)-1
        qp_C=qp_C+1

        abmn(qp_C,1)=ab(qp_I,1)
        abmn(qp_C,2)=ab(qp_I,2)
        abmn(qp_C,3)=ab(qp_I,2)-1
        abmn(qp_C,4)=ab(qp_I,2)+ab_ap
        qp_C=qp_C+1


        do qp_j=1,opt_ch-2

            abmn(qp_C,1)=ab(qp_I,1)
            abmn(qp_C,2)=ab(qp_I,2)
            abmn(qp_C,3)=ab(qp_I,2)+ab_ap*qp_j
            abmn(qp_C,4)=ab(qp_I,2)+ab_ap*qp_j+ab_ap
             if (abmn(qp_C,4)<=channels) then
            qp_C=qp_C+1
            end if

        end do
        end if


      if(ab(qp_I,1)>((opt_ch_q+1-2)*ab_ap)) then

          abmn(qp_C,1)=ab(qp_I,1)
        abmn(qp_C,2)=ab(qp_I,2)
        abmn(qp_C,3)=ab(qp_I,2)-1
        abmn(qp_C,4)=ab(qp_I,1)+1
        qp_C=qp_C+1

        abmn(qp_C,1)=ab(qp_I,1)
        abmn(qp_C,2)=ab(qp_I,2)
        abmn(qp_C,3)=ab(qp_I,1)+1
        abmn(qp_C,4)=ab(qp_I,1)-ab_ap
        qp_C=qp_C+1

        do qp_j=1,opt_ch-2

            abmn(qp_C,1)=ab(qp_I,1)
            abmn(qp_C,2)=ab(qp_I,2)
            abmn(qp_C,3)=ab(qp_I,1)-ab_ap*qp_j
            abmn(qp_C,4)=ab(qp_I,1)-ab_ap*qp_j-ab_ap
             if (abmn(qp_C,4)>0) then
           qp_C=qp_C+1
           end if

        end do
        end if

   end if

   if (ab_ap>k) then  !This represents the case of large apertures. We want only WS and not DD

    kk=max(floor((-opt_ch+(opt_ch**2+4*ab_ap)**0.5)/2),2) !Second order equation to define the optimized scaling for WS MN electrodes

   IF(mod(ab_ap,2)==0 .and. ab_ap>=channels/2) then !A pair of conditions to avoid consider all all the apertures AB possible, rejecting the unnecessary and focusing on investigation depth.

            abmn(qp_C,1)=ab(qp_I,1)
            abmn(qp_C,2)=ab(qp_I,2)
            abmn(qp_C,3)=ab(qp_I,1)+kk
            abmn(qp_C,4)=ab(qp_I,2)-kk
            qp_C=qp_C+1

        do qp_j=2,opt_ch

            abmn(qp_C,1)=ab(qp_I,1)
            abmn(qp_C,2)=ab(qp_I,2)
            abmn(qp_C,3)=abmn(qp_C-1,4)
            if(mod(qp_j,2)==0) then
            abmn(qp_C,4)=abmn(qp_C-1,3)+kk
            else
            abmn(qp_C,4)=abmn(qp_C-1,3)-kk
            end if

            if (abmn(qp_C,4)==abmn(qp_C,3))then
             abmn(qp_C,4)=abmn(qp_C,4)+1
            end if

            if (abmn(qp_C,4)>0 .and. abmn(qp_C,3)>0 .and. abmn(qp_C,4)<=channels .and. abmn(qp_C,3)<=channels) then
           qp_C=qp_C+1
           end if

        end do


    END IF

   end if

    end do

    !Check last quadripole because it could have gone wrong.

    if(abmn(qp_C,3)==abmn(qp_C,4) .or. abmn(qp_C,3)==abmn(qp_C,2) .or. abmn(qp_C,3)==abmn(qp_C,1)) then
    qp_C=qp_C-1
    end if

    if(abmn(qp_C,3)==abmn(qp_C,4) .or. abmn(qp_C,4)==abmn(qp_C,2) .or. abmn(qp_C,4)==abmn(qp_C,1)) then
    qp_C=qp_C-1
    end if

    if(abmn(qp_C,3)<1 .or. abmn(qp_C,4)<1 .or. abmn(qp_C,3)<channels .or. abmn(qp_C,4)<channels) then
    qp_C=qp_C-1
    end if

END SUBROUTINE

SUBROUTINE apertures(channels,ab,max_ap)
    integer, intent (in)  :: channels              ! input
    integer, intent (out) :: max_ap ! output
    integer, dimension(100000,2), intent (out) :: ab ! output
    INTEGER::ap_I,ap_C,ap_J,ap_O

    ap_C=1

    DO ap_I=1,(channels-1)

    ap_J=1
    ap_O=1
    DO WHILE(ap_O==1)

    ab(ap_C,1)=ap_J
    ab(ap_C,2)=ap_J+ap_I

    IF (ab(ap_C,2)==channels) THEN
    ap_O=0
    END IF

    ap_C=ap_C+1
    ap_J=ap_J+1

    END DO

    END DO
max_ap=ap_C-1

END SUBROUTINE

