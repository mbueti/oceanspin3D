#!/bin/sh

PS4='+ gfdl_pre_ocean_sortvit.sh line $LINENO: '

# This script is used to filter through the TC vitals archive.

#eliminate 4-digit years
awk '
{
  yycheck = substr($0,20,2)
  if ((yycheck == 20 || yycheck == 19) && (length($4) == 8)) {
    printf ("%s%s\n",substr($0,1,19),substr($0,22))
  }
  else {
    printf ("%s\n",$0)
  }
} ' | \
#sort by date and hour
sort -n +3 -4 +4 -5 | \
#eliminate repeating dates
awk '
{
  line[NR]=$0
}
END {
  n=0
  for (i=1; i<=NR; i++) {
    j=i
    while (substr(line[j],20,11)==substr(line[j+1],20,11))
      j++
      cmp = length(line[i])
      # Among lines with the same date, select the longest.  Also, by 
      # using a ">=" in the if statement below (as opposed to a ">")
      # we will be taking the LAST, correct, full-length tcvitals 
      # record for this time & storm.  This is important, because the
      # last one issued is always the correct one and invalidates all
      # previous tcvitals for that time & storm (this would be in
      # case a TPC forecaster accidentally messes up and sends out a 
      # tcvitals with incorrect/old info and then sends out a second
      # tcvitals as a replacement).
      ii = i
      for (i1=i+1; i1<=j; i1++)
        if (length(line[i1]) >= cmp) { cmp=length(line[i1]); ii=i1 }
        n++
        lineout[n]=line[ii]
        i=j
  }
  #reformat
  for (i=1; i<=n; i++)
    printf("%s%s%s%s%s\n",substr(lineout[i],1,19),substr(lineout[i],20,15),\
    substr(lineout[i],36,5),substr(lineout[i],42,50),substr(lineout[i],94,20))
} ' | \
#remove data inconsistencies
awk '
{
  n=split($0,field)
  for (i=15; i<=n; i++) {
    if ( field[i] < field[14] ) field[i]="-999"
    if ( field[i] == 999 )      field[i]="-999"
  }
  printf("%s",substr($0,1,69))
  for (i=15; i<=n; i++)
    printf(" %s",field[i])
  printf("\n")
} '

