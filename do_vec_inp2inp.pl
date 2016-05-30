# Comments:

# convert inp files to vector (extract the coefficient values at the nodes)



$N1 = $ARGV[0];
$N2 = $ARGV[1];

print $N1; print "\n";
print $N2; print "\n";


for ($i == $N1; $i <= $N2; $i++)
{
	
	$inputfile = sprintf('%s%d%s', 'object',$i,'_trunc.m');
#	$inputfile2 = sprintf('%s%d%s', '/home/thanhnguyen/UNCC_2012/programs/final_results_parameter_files/grid_obj',$i,'.inp');
	$inputfile2 = sprintf('%s%d%s', 'grid_obj',$i,'.inp');
	$outputfile = sprintf('%s%d%s', 'object',$i,'_final.inp');
	print $inputfile; print ", ";	print $inputfile2; print "\n";
  
	`   /home/tnguy152/UNCC_2012/programs/GCA/vec_inp2inp $inputfile $inputfile2 $outputfile `;


} 

