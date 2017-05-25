$!  Decide whether or not this is a VAX or an Alpha
$!
$ IF "''F$SEARCH("SYS$SYSTEM:VAXVMSSYS.PAR")'" .EQS. ""
$ THEN
$       OPTION = "EXTERNAL.OPT_ALPHA"
$ ELSE
$       OPTION = "EXTERNAL.OPT_VAX"
$ ENDIF
$!
$ VERIFY = F$VERIFY(1)
$ CC REARRANGE /DEFINE="INCLUDE=""REARRANGE.H_GENERIC"""
$ FORTRAN /EXTEND_SOURCE FMEDIAN.F
$ LINK/SHARE=EXTERNAL.EXE REARRANGE,FMEDIAN,'OPTION'/OPT
$ DELETE *.OBJ;*
$ PURGE EXTERNAL.EXE
$ RENAME EXTERNAL.EXE ;1
$ VERIFY = F$VERIFY(VERIFY)