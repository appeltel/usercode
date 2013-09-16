#!/usr/bin/perl


open( TMPL,"ana_PBS_template") or die $!;
@pbs = <TMPL>;
close TMPL;

$seq = 0;
for( $seq=1; $seq <= 50; $seq++)
{
   open( PBS, ">ana" . "_$seq" . ".pbs");

   foreach $line (@pbs)
   {
     $t = $line;
     $t =~ s/#_SEQ_#/$seq/g;
     $f_seq = 50*$seq + 1;
     $t =~ s/#_SEQ50_#/$f_seq/g;
     print PBS $t;
   }
   close PBS;
}

