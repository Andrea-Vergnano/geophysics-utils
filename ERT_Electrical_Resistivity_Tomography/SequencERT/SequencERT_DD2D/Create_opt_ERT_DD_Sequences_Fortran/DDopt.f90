! This program creates a dipole dipole sequence, asking the user how many channels and how many roll-along channels he wish to employ,
!and how many measurements can be acquired at the same time by the instruments (with a logic that suits the IRIS Syscal Pro instrument,
!which can measure up to 10 channels at the same time if the measurements are like: ABMN: 1 2 3 4; 1 2 4 5; 1 2 5 10; 1 2 10 11...)

program DD_sequence

implicit none


INTEGER :: ap, i, j, meas_count, channels,max_ap,min_opt_ch,opt_ch,roll_count,is_present,roll_ch,groupcount,indmin,indmax
  integer, dimension(100000) ::A,B,M,N,Ar,Br,Mr,Nr,group,groupmin,Aord,Bord,Mord,Nord
  integer :: io, ioEl, groupminA,groupminB,groupminM,groupminN
  integer ::ioroll,iorollEl
character(10)::str
meas_count=0
roll_count=0
channels=48
max_ap=10
opt_ch=10

roll_ch=24

! ask the user for input values:


print *, 'This program produces an optimized DD sequence for ERT surveys:'

print *, 'How many channels you have?:'
  read(*,*) channels

  print *, 'How many channels you roll along?:'
  read(*,*) roll_ch

  print *, 'How many channels does your instrument measure at the same time?:'
  read(*,*) opt_ch


  min_opt_ch=max(opt_ch-3,1)


! Build the standard Dipole dipole sequence

DO ap=1,max_ap

if( (channels-min_opt_ch*ap-ap*2)>0 ) then

  do i = 1,(channels-min_opt_ch*ap-ap*2)

 do j = 1,opt_ch


     if(i<channels .and. i+ap < channels .and. i+ap+j*ap < channels .and. i+ap+j*ap+ap<channels) then
     meas_count=meas_count+1
      A(meas_count)=i
      B(meas_count)=i+ap
      M(meas_count)=i+ap+j*ap
      N(meas_count)=i+ap+j*ap+ap
      end if


    end do

    end do
end if

END DO


! add the reverse sequence

DO ap=1,max_ap

if( (channels-min_opt_ch*ap-ap*2)>0 ) then

  do i = 1,(channels-min_opt_ch*ap-ap*2)

 do j = 1,opt_ch

          if(channels+1-i>0 .and. channels+1-(i+ap)>0 .and. channels+1-(i+ap+j*ap)>0 .and. channels+1-(i+ap+j*ap+ap)>0) then
        meas_count=meas_count+1
      A(meas_count)=channels+1-i
      B(meas_count)=channels+1-(i+ap)
      M(meas_count)=channels+1-(i+ap+j*ap)
      N(meas_count)=channels+1-(i+ap+j*ap+ap)
        end if

    end do

    end do
end if

END DO


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


!save csv formatted for resipy forward.


open(newunit=io, file="DDopt.csv", status="new", action="write")

write(io,fmt='(G0.1)') 'A,B,M,N'

do i =1,meas_count
write( io,fmt='(*(G0.1,:,",",:))') Aord(i),Bord(i),Mord(i),Nord(i)
end do

close(io)

!save txt formatted for electre pro
open(newunit=ioEl, file="DDopt.txt", status="new", action="write")

write(ioEl,fmt='(G0.1)') '# X Y Z'

do i =1,channels
write( ioEl,fmt='(*(G0.1,:,Tr1,:))') i, i-1,0,0


end do


write(ioEl,fmt='(G0.1)') '# A B M N'

do i =1,meas_count

write( ioEl,fmt='(*(G0.1,:,Tr1,:))') i,Aord(i),Bord(i),Mord(i),Nord(i)

end do

close(ioEl)



!compute the roll along sequence, equal to the normal sequence but without those measurements already performed in the previous non-roll measurement.


A=Aord
B=Bord
M=Mord
N=Nord

outer:do i=1,meas_count
is_present=0
inner:do j=1,meas_count

if(i/=j .and. A(i)+roll_ch==A(j) .and. B(i)+roll_ch==B(j) .and. M(i)+roll_ch==M(j) .and. N(i)+roll_ch==N(j)) then
is_present=1
exit inner
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

!save csv formatted for resipy forward.

open(newunit=ioroll, file="DDopt_roll.csv", status="new", action="write")

write(ioroll,fmt='(G0.1)') 'A,B,M,N'
do i =1,roll_count

write( ioroll,fmt='(*(G0.1,:,",",:))') Ar(i),Br(i),Mr(i),Nr(i)

end do

close(ioroll)


!save txt formatted for electre pro
open(newunit=iorollEl, file="DDopt_roll.txt", status="new", action="write")

write(iorollEl,fmt='(G0.1)') '# X Y Z'

do i =1,channels
write( iorollEl,fmt='(*(G0.1,:,Tr1,:))') i, i-1,0,0


end do


write(iorollEl,fmt='(G0.1)') '# A B M N'

do i =1,roll_count

write( iorollEl,fmt='(*(G0.1,:,Tr1,:))') i,Ar(i),Br(i),Mr(i),Nr(i)

end do

close(iorollEl)



end program DD_sequence
