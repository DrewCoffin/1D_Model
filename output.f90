MODULE OUTPUTS   

CONTAINS

SUBROUTINE IonOutput(sp, s2p, s3p, op, o2p, longitude, day_char, quantity)
double precision ::sp, s2p, s3p, op, o2p
real             ::longitude 
character(len=4) ::quantity, day_char 
     
  open(unit=101, file=''//quantity//'sp'//day_char//'.dat' , status='unknown', position='append')
  open(unit=102, file=''//quantity//'s2p'//day_char//'.dat', status='unknown', position='append')
  open(unit=103, file=''//quantity//'s3p'//day_char//'.dat', status='unknown', position='append')
  open(unit=104, file=''//quantity//'op'//day_char//'.dat' , status='unknown', position='append')
  open(unit=105, file=''//quantity//'o2p'//day_char//'.dat', status='unknown', position='append')
  write(101,*) longitude, sp
  write(102,*) longitude, s2p
  write(103,*) longitude, s3p
  write(104,*) longitude, op
  write(105,*) longitude, o2p
  close(101)
  close(102)
  close(103)
  close(104)
  close(105)

end subroutine IonOutput

SUBROUTINE IonElecOutput(sp, s2p, s3p, op, o2p, elec, longitude, day_char, quantity)
double precision ::sp, s2p, s3p, op, o2p, elec
real             ::longitude 
character(len=4) ::quantity, day_char 
     
  open(unit=101, file=''//quantity//'sp'//day_char//'.dat' , status='unknown', position='append')
  open(unit=102, file=''//quantity//'s2p'//day_char//'.dat', status='unknown', position='append')
  open(unit=103, file=''//quantity//'s3p'//day_char//'.dat', status='unknown', position='append')
  open(unit=104, file=''//quantity//'op'//day_char//'.dat' , status='unknown', position='append')
  open(unit=105, file=''//quantity//'o2p'//day_char//'.dat', status='unknown', position='append')
  open(unit=106, file=''//quantity//'elec'//day_char//'.dat', status='unknown', position='append')
  write(101,*) longitude, sp
  write(102,*) longitude, s2p
  write(103,*) longitude, s3p
  write(104,*) longitude, op
  write(105,*) longitude, o2p
  write(106,*) longitude, elec
  close(101)
  close(102)
  close(103)
  close(104)
  close(105)
  close(106)

end subroutine IonElecOutput

subroutine OtherOutput(val, longitude, day_char, quantity)
real             ::longitude, val
character(len=4) ::quantity, day_char

  open(unit=101, file=''//quantity//day_char//'.dat', status='unknown', position='append')
  write(101,*) longitude, val
  close(101)

end subroutine OtherOutput

END MODULE
