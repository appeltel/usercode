@skiplist = ();

for( $i = 1; $i<51; $i++)
{

  $good = 1;

  foreach $skip (@skiplist)
  {
    if ($i == $skip ) { $good = 0; }
  }

  if ( $good )
  {
    system("qsub ana_$i.pbs");
  }
}
