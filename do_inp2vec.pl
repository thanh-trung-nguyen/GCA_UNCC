# Comments:

# convert inp files to vector (extract the coefficient values at the nodes)



$N1 = $ARGV[0];
$N2 = $ARGV[1];

print $N1; print "\n";
print $N2; print "\n";


for ($i == $N1; $i <= $N2; $i++)
{
	
	$inputfile = sprintf('%s%d%s', 'object',$i,'.inp');
	$outputfile = sprintf('%s%d%s', 'object',$i,'.m');
print $inputfile; print "\n";
  
	`   /home/tnguy152/UNCC_2012/programs/C_lib/inp2vector $inputfile $outputfile `;


} 

