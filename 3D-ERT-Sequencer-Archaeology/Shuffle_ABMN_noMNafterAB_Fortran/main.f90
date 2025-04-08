program Shuffle_ABMN_noMNafterAB
  implicit none
  integer, parameter :: n = 10000
integer,dimension(100000,4)::abmn,abmnO
  integer ii,io,qp_C,channels,ioEl,i

!!! data input must be named data.csv, have 4 columns of integer numbers and no header


    qp_C=-1
  open (17, file='data.csv', status='old')
  read (17, '()') !
  do ii = 1, n
  qp_C=qp_C+1
    read (17, *,iostat=io) abmn(ii,1), abmn(ii,2),abmn(ii,3),abmn(ii,4)
        IF (io/=0) EXIT
    !print*,abmn(ii,1), abmn(ii,2),abmn(ii,3),abmn(ii,4)
  end do

  close (17)

    print*,qp_C


  channels=maxval(abmn)
  print*,channels


    !Make it optimized for not making potential electrodes just after current electrodes
call reorderEleCP(abmn,qp_C,abmnO,channels) !abmnO is the output of this subroutine
abmn(1:qp_C,:)=abmnO(1:qp_C,:)


   !Save the results in a format for ElectrePro

open(newunit=ioEl, file="data_opt.txt", status="new", action="write")

!write(ioEl,fmt='(G0.1)') '# X Y Z'

!do i =1,channels
!write( ioEl,fmt='(*(G0.1,:,Tr1,:))') i, i-1,0,0


!end do


write(ioEl,fmt='(G0.1)') '# A B M N'

do i =1,qp_C-1

write( ioEl,fmt='(*(G0.1,:,Tr1,:))') i,abmn(i,1),abmn(i,2),abmn(i,3),abmn(i,4)

end do

close(ioEl)



end program Shuffle_ABMN_noMNafterAB



SUBROUTINE reorderEleCP(abmn,qp_C,abmnO,channels)
    integer, dimension(100000,4), intent (in) :: abmn ! input
    integer, dimension(100000,4), intent (out) :: abmnO ! output
    integer,intent(in)::qp_C,channels

INTEGER::groupcount,groupminA,groupminB,groupminM,groupminN,indmin,indmax,i,di,meas_count,checkdepth,min_groupdi
INTEGER::GrShcount,ri,iter,k,maxiter,max_correction
real::randshuffle,cost
INTEGER,dimension(100000)::group,groupdi,A,B,M,N,Aord,Bord,Mord,Nord,GrSh
integer, dimension(100000,4)::abmntemp

checkdepth=10
min_groupdi=1
maxiter=3000 !!! to be changed if you see convergence in the values of the cost function printed during run
max_correction=3 !!!!! to be changed according to how many electrodes you have. 3 is ok for 48, 4 for 72...


abmnO=abmn
abmntemp=abmn


meas_count=qp_C



! order the sequence in a way one potential electrodes are not selected among those just been current electrodes


!creates the group variable


do iter=1,maxiter

groupcount=1
group(1)=1

! calculate the first index of each optimized group of measurements

do i=2, meas_count

if(abmnO(i,1)/=abmnO(i-1,1) .or.abmnO(i,2)/=abmnO(i-1,2) .or. abmnO(i,3)/=abmnO(i-1,4)) then   !consider also that 20 measurements could happened with the same AB, not only 10.


groupcount=groupcount+1
group(groupcount)=i
end if
end do

group(groupcount+1)=meas_count+1

! group 1 quindi va da 1 a group(2)-1, e così via; group finale è group(groupcount):meas_count
!di ciascun gruppo devo calcolare di=j-i secondo wilkinson2012.



do i=2,groupcount !check the #checkdepth groups before the current group, if they contain the current electrode.

di=checkdepth
groupdi(1)=checkdepth

check: do di= 1,min(checkdepth,i)

if(any(abmnO(group(i):(group(i+1)-1),3:4)==(abmnO(group(i-di),1)))) then
exit check
end if

if(any(abmnO(group(i):(group(i+1)-1),3:4)==(abmnO(group(i-di),2)))) then
exit check
end if

end do check

!at this point di should be the di of Wilkinson 2012, and it has maximum value of checkdepth, which is fine for us.

groupdi(i)=di

end do



!evaluate cost function which should be minimized after some iterations

cost=0.0
do k=1,groupcount
cost=cost+1.0/(groupdi(k)**2)
end do

!do i=1,groupcount
!print *, group(i),groupdi(i)!,NEW_LINE('a')
!end do

print *,cost, groupcount

!The groups of measurements with least groupdi have to be moved!

!choose at random the least groupdi groups (not the first #checkdepth)

!min_groupdi=MINVAL(groupdi(1:groupcount))


!if there are too few min, we want to consider also the minval+1
flagmin=0
GrShcount=0
do i = 1,groupcount
if(groupdi(i)<=min_groupdi )then ! in case there is only one minimum the algorithm could get stuck
GrShcount=GrShcount+1
end if
end do

!locate the worst groups which have to be reshuffled


if(GrShcount<30 .and. min_groupdi<max_correction)then! checkdepth-1)then! put <10

!!!!!

min_groupdi=min_groupdi+1
end if
!print*,'flag',flagmin
GrShcount=0

!$OMP PARALLEL
do i = 1,groupcount-1
if(groupdi(i)<=min_groupdi)then
!if(groupdi(i)==min_groupdi .or.(groupdi(i)==min_groupdi+1 .and. flagmin==1))then ! in case there is only one minimum the algorithm could get stuck
GrShcount=GrShcount+1
GrSh(GrShcount)=i
end if
end do
!$OMP END PARALLEL
!shuffle those groups according to Fisher and Yates' original method (wikipedia)

do i =1,GrShcount-1

call random_number(randshuffle)
!if(randshuffle<0.01)then
ri=floor(randshuffle*(GrShcount-i))+1+i ! this -i and +1 serves to fetch random numbers only between those not already used

abmnO(group(GrSh(i)):group(GrSh(i)+1)-1,:)=abmntemp(group(GrSh(ri)):group(GrSh(ri)+1)-1,:)
abmnO(group(GrSh(ri)):group(GrSh(ri)+1)-1,:)=abmntemp(group(GrSh(i)):group(GrSh(i)+1)-1,:)

!ri=floor(randshuffle*(groupcount-checkdepth))+1.0+checkdepth ! this -i and +1 serves to fetch random numbers only between those not already used

!abmnO(group(GrSh(i)):group(GrSh(i)+1)-1,:)=abmntemp(group(ri):group(ri+1)-1,:)
!abmnO(group(ri):group(ri+1)-1,:)=abmntemp(group(GrSh(i)):group(GrSh(i)+1)-1,:)


abmntemp=abmnO
!end if

end do

end do

!print *, 'minval', MINVAL(groupdi(1:groupcount))
!print *, 'loc',GrSh(1:GrShcount)

!put the results in Aord, Bord, Mord, Nord arrays

do i=1,groupcount
print *, group(i),groupdi(i)!,NEW_LINE('a')
end do

print*,'min_groupdi',min_groupdi

END SUBROUTINE

